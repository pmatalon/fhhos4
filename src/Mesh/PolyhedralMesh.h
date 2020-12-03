#pragma once
#include <set>
#include <list>
#include <algorithm>
#include "Mesh.h"
#include "2D/TriangularElement.h"
#include "2D/QuadrilateralElement.h"
#include "2D/PolygonalElement.h"
#include "3D/TetrahedralElement.h"
#include "3D/ParallelepipedElement.h"
#include "Agglo.h"
using namespace std;

template <int Dim>
class PolyhedralMesh : public Mesh<Dim>
{
protected:
	vector<QuadrilateralElement> _quadrilateralElements;
	vector<TriangularElement> _triangularElements;
	vector<TetrahedralElement> _tetrahedralElements;
	vector<ParallelepipedElement> _parallelepipedElements;

	vector<Edge> _edgeFaces;
	vector<TriangularFace> _triangularFaces;
private:
	double _h = 0;
	double _regularity = 0;
protected:
	string _geometryDescription;
public:
	PolyhedralMesh() : Mesh<Dim>()
	{}

	PolyhedralMesh(string geometryDescription) : Mesh<Dim>()
	{
		_geometryDescription = geometryDescription;
	}

	virtual string Description() override
	{
		if (Dim == 2)
			return "Polygonal";
		return "Polyhedral";
	}
	virtual string FileNamePart() override
	{
		return "polyhedral";
	}
	virtual string GeometryDescription() override
	{
		return _geometryDescription;
	}

	virtual double H() override
	{
		if (_h == 0)
		{
			for (Element<Dim>* e : this->Elements)
				_h = max(_h, e->Diameter());
		}
		return _h;
	}

	virtual double Regularity() override
	{
		if (_regularity == 0)
		{
			for (Element<Dim>* e : this->Elements)
				_regularity = min(_regularity, e->Regularity());
		}
		return _regularity;
	}

	virtual size_t MemoryUsage() override
	{
		size_t verticesUsage = this->Vertices.size() * sizeof(Vertex);
		size_t verticesPointers = this->Vertices.size() * sizeof(Vertex*);

		size_t elementsUsage = _quadrilateralElements.size() * sizeof(QuadrilateralElement)
			+ _triangularElements.size() * sizeof(TriangularElement)
			+ _tetrahedralElements.size() * sizeof(TetrahedralElement)
			+ _parallelepipedElements.size() * sizeof(ParallelepipedElement);
		if (Dim == 2)
			elementsUsage += this->Elements.size() * 3 * sizeof(Face<Dim>*); // at least 3 faces
		else if (Dim == 2)
			elementsUsage += this->Elements.size() * 4 * sizeof(Face<Dim>*); // at least 4 faces
		size_t elementsPointers = this->Elements.size() * sizeof(Element<Dim>*);

		size_t oneFaceUsage = 0;
		if (Dim == 2)
			oneFaceUsage = sizeof(Edge);
		else if (Dim == 2)
			oneFaceUsage = sizeof(TriangleIn3D);
		size_t facesUsage = this->Faces.size() * oneFaceUsage;
		size_t facesPointers = this->Faces.size() * sizeof(Face<Dim>*);

		return verticesUsage + verticesPointers + elementsUsage + elementsPointers + facesUsage + facesPointers;
	}

	virtual void CoarsenMesh(CoarseningStrategy strategy, int coarseningFactor) override
	{
		if (this->CoarseMesh)
			return;
		
		if (strategy == CoarseningStrategy::AgglomerationCoarsening)
			CoarsenByAgglomerationByVertexRemoval();
		else if (strategy == CoarseningStrategy::AgglomerationCoarseningByMostCoplanarFaces)
			CoarsenByAgglomerationByMostCoplanarFaces();
		else if (strategy == CoarseningStrategy::AgglomerationCoarseningByClosestCenter 
			  || strategy == CoarseningStrategy::AgglomerationCoarseningByClosestFace 
			  || strategy == CoarseningStrategy::AgglomerationCoarseningByLargestInterface)
			CoarsenByAgglomerationByPairs(strategy);
		else if (strategy == CoarseningStrategy::AgglomerationCoarseningBySeedPoints)
			CoarsenByAgglomerationBySeedPoints();
		else if (strategy == CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours)
			CoarsenByAgglomerationByFaceNeighbours();
		else if (strategy == CoarseningStrategy::AgglomerationCoarseningByVertexNeighbours)
			CoarsenByAgglomerationByVertexNeighbours();
		else if (strategy == CoarseningStrategy::FaceCoarsening)
			FaceCoarsening();
		else
			Mesh<Dim>::CoarsenMesh(strategy, coarseningFactor);
	}

protected:
	virtual void FinalizeCoarsening() override
	{
		// Compute triangulation and bounding box of the polygons
		ElementParallelLoop<Dim> parallelLoop(this->CoarseMesh->Elements);
		parallelLoop.Execute([](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				// TODO: Make Init() available for all element classes
				PolygonalElement* p = dynamic_cast<PolygonalElement*>(e);
				if (p)
				{
					p->ComputeMinimalTriangulation();
					p->ComputeBoundingBox();
				}
			});

		Mesh<Dim>::FinalizeCoarsening();
	}

private:
	//-------------------------------------------------//
	//                                                 //
	//        Coarsening by most coplanar faces        //
	//                                                 //
	//-------------------------------------------------//

	void CoarsenByAgglomerationByMostCoplanarFaces()
	{
		PolyhedralMesh<Dim>* coarseMesh = new PolyhedralMesh<Dim>();
		this->InitializeCoarsening(coarseMesh);
		coarseMesh->ComesFrom.CS = CoarseningStrategy::AgglomerationCoarseningByMostCoplanarFaces;

		//BigNumber elementToAnalyze = 999999999999;

		// Associate to each vertex the list of faces it connects
		map<Vertex*, vector<Face<Dim>*>> vertexFaceMap = this->BuildVertexFaceMap();

		for (Vertex* vertex : this->Vertices)
		{
			vector<Face<Dim>*> vertexFaces = vertexFaceMap[vertex];
			if (vertexFaces.empty()) // this vertex was eliminated because a neighbouring one has been processed
				continue;

			if (vertexFaces.size() < 3)
			{
				vertexFaceMap[vertex].clear();
				continue;
			}
			
			// Process only if interior vertex
			/*bool isOnBoundary = false;
			for (Face<Dim>* f : vertexFaces)
			{
				if (f->IsDomainBoundary)
				{
					isOnBoundary = true;
					break;
				}
			}
			if (isOnBoundary)
				continue;*/

			//--------------------------------------------------------------------------------------------------//
			// Find the couple of faces (f1,f2) that has the largest angle (in [Pi/2; 3Pi/2])                   //
			// It is identified by the maximum, in absolute value, of the inner products of the normal vectors. //
			//--------------------------------------------------------------------------------------------------//

			double maxNormalInnerProd = 0;
			Face<Dim>* f1 = nullptr;
			Face<Dim>* f2 = nullptr;
			for (int i = 0; i < vertexFaces.size(); i++)
			{
				Face<Dim>* fi = vertexFaces[i];
				if (fi->IsRemovedOnCoarserGrid)
					continue;

				DimVector<Dim> ni = fi->Element1->OuterNormalVector(fi);
				for (int j = i + 1; j < vertexFaces.size(); j++)
				{
					Face<Dim>* fj = vertexFaces[j];
					if (fj->IsRemovedOnCoarserGrid)
						continue;
					DimVector<Dim> nj = fj->Element1->OuterNormalVector(fj);

					// Verify that the angle is obtuse
					if (Vect<Dim>(vertex, fi->Center()).dot(Vect<Dim>(vertex, fj->Center())) < 0)
					{
						double normalInnerProd = abs(ni.dot(nj));
						if (normalInnerProd > maxNormalInnerProd)
						{
							f1 = fi;
							f2 = fj;
							maxNormalInnerProd = normalInnerProd;
						}
					}
				}
			}

			if (!f1 && !f2) // No faces can be merges because there is no obtuse angle
				continue;

			//------------------------------------------------------------------------------------------------------//
			// Element agglomeration:                                                                               //
			// on each side of the macro-face f1-f2, the elements connected to the current vertex are agglomerated. //
			//------------------------------------------------------------------------------------------------------//
			
			vector<Element<Dim>*> agglomerate1;
			Element<Dim>* e = f1->Element1;
			agglomerate1.push_back(e);
			while (!e->HasFace(f2))
			{
				for (Face<Dim>* f : e->Faces)
				{
					if (f != f1 && f->IsIn(vertexFaces))
					{
						Element<Dim>* neighbour = f->GetNeighbour(e);
						if (!neighbour->IsIn(agglomerate1))
						{
							agglomerate1.push_back(neighbour);
							e = neighbour;
							break;
						}
					}
				}
			}

			Element<Dim>* coarseElement1 = coarseMesh->AgglomerateFineElements(agglomerate1);

			// Do not process the neighbouring vertices
			for (Element<Dim>* ea : agglomerate1)
			{
				for (auto v : ea->Vertices())
				{
					if (v != vertex)
						vertexFaceMap[v].clear();
				}
			}

			if (!f1->IsDomainBoundary)
			{
				vector<Element<Dim>*> agglomerate2;
				e = f1->Element2;
				agglomerate2.push_back(e);
				while (!e->HasFace(f2))
				{
					for (Face<Dim>* f : e->Faces)
					{
						if (f != f1 && f->IsIn(vertexFaces))
						{
							Element<Dim>* neighbour = f->GetNeighbour(e);
							if (!neighbour->IsIn(agglomerate2))
							{
								agglomerate2.push_back(neighbour);
								e = neighbour;
								break;
							}
						}
					}
				}
				Element<Dim>* coarseElement2 = coarseMesh->AgglomerateFineElements(agglomerate2);

				// Do not process the neighbouring vertices
				for (Element<Dim>* ea : agglomerate2)
				{
					for (auto v : ea->Vertices())
					{
						if (v != vertex)
							vertexFaceMap[v].clear();
					}
				}
			}

			//--------------------------------------------//
			//           Agglomerate f1 and f2            //
			//--------------------------------------------//

			// If f1 and f2 are not totally coplanar, then the nestedness is lost.
			// However, if they are close to coplanar, it's still okay because the coarse mesh can be seen
			// as a slightly perturbed nested mesh.

			FaceCollapsingStatus status = coarseMesh->AgglomerateFineFaces(f1, f2);
			assert(status == FaceCollapsingStatus::Ok);

			coarseMesh->TryCollapseInterfacesMadeOfMultipleFaces(coarseElement1);
		}

		vertexFaceMap.clear();

		//--------------------------------------------------------------------------------------//
		// At this point, there might still be some fine elements that have not been coarsened. //
		// We agglomerate them with their smallest neighbour.                                   //
		//--------------------------------------------------------------------------------------//
		for (Element<Dim>* e : this->Elements)
		{
			if (e->CoarserElement)
				continue;

			Element<Dim>* smallestNeighbour = nullptr;
			double smallestNeighbourSize = 0;
			bool smallestNeighbourIsFine = false;
			for (Face<Dim>* f : e->Faces)
			{
				if (f->IsDomainBoundary)
					continue;
				Element<Dim>* neighbour = f->GetNeighbour(e);
				if (!neighbour->CoarserElement)
				{
					smallestNeighbour = neighbour;
					smallestNeighbourIsFine = true;
					break;
				}
				else
				{
					Element<Dim>* coarseNeigbour = neighbour->CoarserElement;
					double coarseNeigbourSize = coarseNeigbour->Measure();
					if (!smallestNeighbour || coarseNeigbourSize < smallestNeighbourSize)
					{
						smallestNeighbour = coarseNeigbour;
						smallestNeighbourSize = coarseNeigbourSize;
					}
				}
			}

			if (smallestNeighbourIsFine)
				coarseMesh->AgglomerateFineElements({ e, smallestNeighbour });
			else
				coarseMesh->AgglomerateFineElementToCoarse(e, smallestNeighbour);
		}


		//------------------------------------------------------------------------------//
		// We merge the faces which are part of the same interface between two elements //
		//------------------------------------------------------------------------------//
		for (Element<Dim>* ce : coarseMesh->Elements)
			coarseMesh->TryCollapseInterfacesMadeOfMultipleFaces(ce);

		// Set the coarse vertex list
		set<Vertex*> vertices;
		for (auto f : coarseMesh->Faces)
		{
			for (auto v : f->Vertices())
				vertices.insert(v);
		}
		coarseMesh->Vertices = vector<Vertex*>(vertices.begin(), vertices.end());

		this->FinalizeCoarsening();
	}

	//-------------------------------------------------------------------------------------------------------------------------//
	// Browse the neighbours: if the interface with one of them is made of multiple faces, they're collapsed into a single one //
	//-------------------------------------------------------------------------------------------------------------------------//
	void TryCollapseInterfacesMadeOfMultipleFaces(Element<Dim>* e, bool makeItThreadSafe = false)
	{
		if (makeItThreadSafe && !e->Mutex.try_lock())
			return;

		bool retryIfElementChanges = false;
			
		for (Element<Dim>* n : e->Neighbours(false))
		{
			if (makeItThreadSafe && !n->Mutex.try_lock())
				continue;

			FaceCollapsingStatus status = TryCollapseInterfaceBetween(e, n);
			if (status == FaceCollapsingStatus::Ok && retryIfElementChanges)
			{
				if (makeItThreadSafe)
				{
					e->Mutex.unlock();
					n->Mutex.unlock();
				}
				TryCollapseInterfacesMadeOfMultipleFaces(e, makeItThreadSafe);
				return;
			}
			else if (status == FaceCollapsingStatus::ElementFullDegeneration || status == FaceCollapsingStatus::OneElementEmbeddedInConvexHullOfTheOther)
			{
				// We can't collpase the faces, otherwise one element will degenerate. Agglomerating the elements instead.
				Element<Dim>* newE = Agglomerate(e, n);
				if (makeItThreadSafe)
				{
					e->Mutex.unlock();
					n->Mutex.unlock();
				}
				TryCollapseInterfacesMadeOfMultipleFaces(newE, makeItThreadSafe);
				return;
			}
			else if (status == FaceCollapsingStatus::CrossedPolygon || status == FaceCollapsingStatus::ElementPartialDegeneration)
			{
				retryIfElementChanges = true;
			}
			// Any other failing status, abort the collapsing.

			if (makeItThreadSafe)
				n->Mutex.unlock();
		}

		if (e->IsOnBoundary())
			TryCollapseBoundaryFaces(e);

		if (makeItThreadSafe)
			e->Mutex.unlock();
	}

	FaceCollapsingStatus TryCollapseInterfaceBetween(Element<Dim>* e1, Element<Dim>* e2)
	{
		vector<Face<Dim>*> interfaceFaces = e1->InterfaceWith(e2);
		assert(!interfaceFaces.empty());

		if (interfaceFaces.size() == 1)
			return FaceCollapsingStatus::NotEnoughFaces;

		if (e1->IsInSamePhysicalPartAs(e2))
			return TryCollapse(interfaceFaces);
		else
		{
			if (Utils::ProgramArgs.Solver.MG.BoundaryFaceCollapsing == FaceCollapsing::Disabled)
				return FaceCollapsingStatus::NotEnoughFaces;
			else if (Utils::ProgramArgs.Solver.MG.BoundaryFaceCollapsing == FaceCollapsing::OnlyCollinear)
			{
				Interface<Dim> interf(interfaceFaces);
				return TryCollapseCoplanarFaces(interf);
			}
			else if (Utils::ProgramArgs.Solver.MG.BoundaryFaceCollapsing == FaceCollapsing::ByPairs)
				return TryCollapseByPairs(interfaceFaces);
			else if (Utils::ProgramArgs.Solver.MG.BoundaryFaceCollapsing == FaceCollapsing::Max)
				return TryCollapse(interfaceFaces);
			else
				Utils::FatalError("Unmanaged FaceCollapsing mode");
		}
	}

	FaceCollapsingStatus TryCollapseBoundaryFaces(Element<Dim>* e)
	{
		vector<Face<Dim>*> boundaryFaces = e->BoundaryFaces();
		if (boundaryFaces.size() < 2)
			return FaceCollapsingStatus::NotEnoughFaces;

		if (Utils::ProgramArgs.Solver.MG.BoundaryFaceCollapsing == FaceCollapsing::Disabled)
			return FaceCollapsingStatus::NotEnoughFaces;
		else if (Utils::ProgramArgs.Solver.MG.BoundaryFaceCollapsing == FaceCollapsing::OnlyCollinear)
		{
			Interface<Dim> interf(boundaryFaces);
			return TryCollapseCoplanarFaces(interf);
		}
		else if (Utils::ProgramArgs.Solver.MG.BoundaryFaceCollapsing == FaceCollapsing::ByPairs)
			return TryCollapseByPairs(boundaryFaces);
		else if (Utils::ProgramArgs.Solver.MG.BoundaryFaceCollapsing == FaceCollapsing::Max)
			return TryCollapse(boundaryFaces);
		else
			Utils::FatalError("Unmanaged FaceCollapsing mode");
	}

	FaceCollapsingStatus TryCollapseCoplanarFaces(const Interface<Dim>& interf)
	{
		list<set<Face<Dim>*>> coplanarSubsets = interf.CoplanarSubsets();
		if (coplanarSubsets.empty())
			return FaceCollapsingStatus::NotEnoughFaces;
		for (set<Face<Dim>*> subset : coplanarSubsets)
		{
			Interface<Dim> subInterf(vector<Face<Dim>*>(subset.begin(), subset.end()));
			Collapse(subInterf);
		}
		return FaceCollapsingStatus::Ok;
	}



	map<Vertex*, vector<Face<Dim>*>> BuildVertexFaceMap()
	{
		// Init the map with an empty list of faces for each vertex
		map<Vertex*, vector<Face<Dim>*>> vertexFaces;
		for (Vertex* v : this->Vertices)
			vertexFaces.insert({ v, vector<Face<Dim>*>() });

		// Add the faces to the list of their respective vertices
		for (Face<Dim>* f : this->Faces)
		{
			for (Vertex* v : f->Vertices())
				vertexFaces[v].push_back(f);
		}

		// TODO: sort the lists in direct order
		return vertexFaces;
	}

	//-----------------------------------------//
	//                                         //
	//        Coarsening by seed points        //
	//                                         //
	//-----------------------------------------//

	map<Element<Dim>*, Element<Dim>*> _seedElements;
	Element<Dim>* GetSeed(Element<Dim>* e)
	{
		return _seedElements[e];
	}
	void AddSeedAssociation(Element<Dim>* e, Element<Dim>* seed)
	{
		_seedElements[e] = seed;
	}
	Element<Dim>* ClosestSeed(Element<Dim>* e)
	{
		Element<Dim>* closestSeed = nullptr;
		double closestSeedDistance = -1;

		for (Face<Dim>* f : e->Faces)
		{
			if (f->IsDomainBoundary)
				continue;

			Element<Dim>* neighbour = f->GetNeighbour(e);
			if (neighbour->CoarserElement)
			{
				Element<Dim>* neighbourSeed = GetSeed(neighbour);
				double distance = Vect<2>(e->Center(), neighbourSeed->Center()).norm();
				if (!closestSeed || distance < closestSeedDistance)
				{
					closestSeed = neighbourSeed;
					closestSeedDistance = distance;
				}
			}
		}

		return closestSeed;
	}

	void CoarsenByAgglomerationBySeedPoints()
	{
		PolyhedralMesh<Dim>* coarseMesh = new PolyhedralMesh<Dim>();
		this->InitializeCoarsening(coarseMesh);
		coarseMesh->ComesFrom.CS = CoarseningStrategy::AgglomerationCoarseningBySeedPoints;

		for (Element<Dim>* e : this->Elements)
		{
			if (CanBeSeed(e))
				coarseMesh->AgglomerateNeighbours(e);
		}

		// For the remaining elements, they're aggregated with the closest seed
		for (int i = 0; i < 2; i++)
		{
			for (Element<Dim>* e : this->Elements)
			{
				if (e->CoarserElement)
					continue;

				Element<Dim>* seed = ClosestSeed(e);
				if (seed)
				{
					coarseMesh->AgglomerateFineElementToCoarse(e, seed->CoarserElement);
					AddSeedAssociation(e, seed);
				}
			}
		}

		_seedElements.clear();
		this->FinalizeCoarsening();
	}

	bool CanBeSeed(Element<Dim>* e)
	{
		if (e->CoarserElement)
			return false;

		for (Face<Dim>* f : e->Faces)
		{
			if (f->IsDomainBoundary)
				return false;

			Element<Dim>* neighbour = f->GetNeighbour(e);
			if (neighbour->CoarserElement)
				return false;
		}
		return true;
	}

	Element<Dim>* AgglomerateNeighbours(Element<Dim>* e)
	{
		Element<Dim>* macroElement = nullptr;
		for (Face<Dim>* f : e->Faces)
		{
			if (f->IsDomainBoundary)
				continue;

			Element<Dim>* neighbour = f->GetNeighbour(e);
			if (neighbour->CoarserElement) // because the interface between e and neighbour can be composed of multiple faces
				continue;

			if (!macroElement)
				macroElement = AgglomerateFineElements(e, neighbour);
			else
				macroElement = AgglomerateFineElementToCoarse(neighbour, macroElement);
			dynamic_cast<PolyhedralMesh<Dim>*>(this->FineMesh)->AddSeedAssociation(neighbour, e);
		}
		dynamic_cast<PolyhedralMesh<Dim>*>(this->FineMesh)->AddSeedAssociation(e, e);
		return macroElement;
	}

	//---------------------------------------------//
	//                                             //
	//        Coarsening by face neighbours        //
	//                                             //
	//---------------------------------------------//

	void CoarsenByAgglomerationByFaceNeighbours()
	{
		PolyhedralMesh<Dim>* coarseMesh = new PolyhedralMesh<Dim>();
		this->InitializeCoarsening(coarseMesh);
		coarseMesh->ComesFrom.CS = CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours;

		// Element agglomeration, 1st pass
		// Agglomerate fine elements together

		bool cancelCoarsening = false;

		// Start by the elements at re-entrant corners
		if (Utils::ProgramArgs.Solver.MG.ReEntrantCornerManagement == ReEntrantCornerMgmt::AgglomerateFirst)
		{
			for (Element<Dim>* e : this->Elements)
			{
				if (!e->HasReEntrantCorner || e->CoarserElement)
					continue;

				vector<Element<Dim>*> availableNeighbours = LockAvailableFaceNeighbours(e, true);
				if (availableNeighbours.empty())
				{
					e->Mutex.unlock();
					/*Utils::Warning("No more agglomeration possible in this mesh or this physical region. Coarsening aborted.");
					cancelCoarsening = true;
					coarseMesh->ExportToMatlab2("/mnt/c/Users/pierr/Desktop/mesh_in_construct.m");
					e->ExportToMatlab("k");
					return;*/
					continue;
				}

				// Agglomeration
				availableNeighbours.push_back(e);
				Element<Dim>* coarseElement = coarseMesh->AgglomerateFineElements(availableNeighbours, true);
				if (!coarseElement)
				{
					Utils::Warning("No more agglomeration possible in this mesh or this physical region. Coarsening aborted.");
					cancelCoarsening = true;
					//return;
				}
				for (Element<Dim>* neighbour : availableNeighbours)
					neighbour->Mutex.unlock();
				if (cancelCoarsening)
					return;
			}
		}


		vector<Element<Dim>*> remainingFineElements = this->Elements;
		random_shuffle(remainingFineElements.begin(), remainingFineElements.end());
		//cout << "\t" << remainingFineElements.size() << " fine elements to coarsen" << endl;

		bool elementsAreAgglomerated = true;
		while (elementsAreAgglomerated)
		{
			elementsAreAgglomerated = false;
			ElementParallelLoop<Dim> parallelLoop(remainingFineElements);
			parallelLoop.Execute([this, coarseMesh, &elementsAreAgglomerated](Element<Dim>* currentElem)
				{
					if (currentElem->CoarserElement)
						return;

					if (!currentElem->Mutex.try_lock())
						return;

					// Lock neighbours
					vector<Element<Dim>*> availableNeighbours = LockAvailableFaceNeighbours(currentElem);
					if (availableNeighbours.empty())
					{
						currentElem->Mutex.unlock();
						return;
					}

					// Agglomeration
					availableNeighbours.push_back(currentElem);
					Element<Dim>* coarseElement = coarseMesh->AgglomerateFineElements(availableNeighbours);
					if (coarseElement)
						elementsAreAgglomerated = true;

					for (Element<Dim>* neighbour : availableNeighbours)
						neighbour->Mutex.unlock();
				});

			if (elementsAreAgglomerated)
			{
				RemoveCoarsenedElements(remainingFineElements);
				//cout << "\t" << remainingFineElements.size() << " fine elements to coarsen" << endl;
			}
		}

		// Element agglomeration, 2nd pass
		// Agglomerate the remaining fine elements with their closest coarse neighbour

		//cout << "\t" << remainingFineElements.size() << " to agglomerate with coarse elements" << endl;

		while (!remainingFineElements.empty())
		{
			ElementParallelLoop<Dim> parallelLoop(remainingFineElements);
			parallelLoop.Execute([this, coarseMesh, &cancelCoarsening](Element<Dim>* currentElem)
				{
					if (cancelCoarsening)
						return;
					assert(!currentElem->CoarserElement);
					Element<Dim>* coarseNeighbourForAggreg = this->FittestCoarseNeighbour(currentElem, CoarseningStrategy::AgglomerationCoarseningByClosestCenter);
					if (!coarseNeighbourForAggreg)
					{
						if (!cancelCoarsening)
						{
							cancelCoarsening = true;
							//coarseMesh->ExportToMatlab2("/mnt/c/Users/pierr/Desktop/mesh_in_construct.m");
							//currentElem->ExportToMatlab("k");
							Utils::Warning("No more agglomeration possible in this mesh or this physical region. Coarsening aborted.");
						}
						return;
					}

					if (coarseNeighbourForAggreg->Mutex.try_lock())
					{
						if (!coarseNeighbourForAggreg->IsDeleted)
						{
							assert(!coarseNeighbourForAggreg->IsDeleted);
							Element<Dim>* coarseElement = coarseMesh->AgglomerateFineElementToCoarse(currentElem, coarseNeighbourForAggreg);
						}
						coarseNeighbourForAggreg->Mutex.unlock();
					}
				});

			if (cancelCoarsening)
			{
				delete coarseMesh;
				this->CoarseMesh = nullptr;
				return;
			}

			RemoveCoarsenedElements(remainingFineElements);
			//cout << "\t" << remainingFineElements.size() << " fine elements to coarsen" << endl;
		}

		// Face collapsing
		ElementParallelLoop<Dim> parallelLoopCollapseFaces(coarseMesh->Elements);
		parallelLoopCollapseFaces.Execute([coarseMesh](Element<Dim>* coarseElement)
			{
				if (!coarseElement->IsDeleted)
					coarseMesh->TryCollapseInterfacesMadeOfMultipleFaces(coarseElement, true);
			});

		for (Face<Dim>* f : coarseMesh->Faces)
		{
			if (!f->IsDeleted && !f->IsDomainBoundary && !f->FinerFaces[0]->HasBeenCoarsened())
				coarseMesh->TryCollapseInterfacesMadeOfMultipleFaces(f->Element1, false);
		}
		
		/*list<Face<Dim>*> uncoarsenedFaces;
		for (Face<Dim>* f : this->Faces)
		{
			if (!f->IsDeleted && !f->IsRemovedOnCoarserGrid && !f->HasBeenCoarsened())
				uncoarsenedFaces.push_back(f);
		}
		cout << (uncoarsenedFaces.size() * 100) / this->Faces.size() << "% unchanged faces." << endl;

		int nFacesAtPhysicalBoundary = 0;
		int nUncorsenedFacesAtPhysicalBoundary = 0;
		for (Face<Dim>* f : this->Faces)
		{
			if (f->IsDeleted || f->IsDomainBoundary)
				continue;
			bool isAtPhysicalBoundary = !f->Element1->IsInSamePhysicalPartAs(f->Element2);
			if (!isAtPhysicalBoundary)
				continue;
			nFacesAtPhysicalBoundary++;
			if (!f->IsRemovedOnCoarserGrid && !f->HasBeenCoarsened())
				nUncorsenedFacesAtPhysicalBoundary++;
		}
		if (nFacesAtPhysicalBoundary > 0)
			cout << (nUncorsenedFacesAtPhysicalBoundary * 100) / nFacesAtPhysicalBoundary << "% unchanged faces at physical boundaries." << endl;*/

		this->FinalizeCoarsening();
	}

	//-----------------------------------------------//
	//                                               //
	//        Coarsening by vertex neighbours        //
	//                                               //
	//-----------------------------------------------//

	void CoarsenByAgglomerationByVertexNeighbours()
	{
		PolyhedralMesh<Dim>* coarseMesh = new PolyhedralMesh<Dim>();
		this->InitializeCoarsening(coarseMesh);
		coarseMesh->ComesFrom.CS = CoarseningStrategy::AgglomerationCoarseningByVertexNeighbours;

		// Associate to each vertex the list of elements it connects
		map<Vertex*, vector<Element<Dim>*>> vertexElements = BuildVertexElementMap();

		for (Element<Dim>* currentElem : this->Elements)
		{
			if (currentElem->CoarserElement)
				continue;

			vector<Element<Dim>*> availableNeighbours = AvailableVertexNeighbours(currentElem, vertexElements);

			if (availableNeighbours.size() > 0)
			{
				availableNeighbours.push_back(currentElem);
				Element<Dim>* coarseElement = coarseMesh->AgglomerateFineElements(availableNeighbours);
				coarseMesh->TryCollapseInterfacesMadeOfMultipleFaces(coarseElement);
			}
			else
			{
				Element<Dim>* coarseNeighbourForAggreg = this->FittestCoarseNeighbour(currentElem, CoarseningStrategy::AgglomerationCoarseningByClosestCenter);
				if (!coarseNeighbourForAggreg)
					Utils::FatalError("Element cannot be aggregated. Weird...");

				Element<Dim>* macroElement = coarseMesh->AgglomerateFineElementToCoarse(currentElem, coarseNeighbourForAggreg);
				coarseMesh->TryCollapseInterfacesMadeOfMultipleFaces(macroElement);
			}

		}

		vertexElements.clear();
		this->FinalizeCoarsening();
	}

	//--------------------------------------------//
	//                                            //
	//        Coarsening by vertex removal        //
	//                                            //
	//--------------------------------------------//

	//-----------------//
	// NOT FINISHED!!! //
	//-----------------//
	void CoarsenByAgglomerationByVertexRemoval()
	{
		PolyhedralMesh<Dim>* coarseMesh = new PolyhedralMesh<Dim>();
		this->InitializeCoarsening(coarseMesh);

		// Associate to each vertex the list of elements it connects
		map<Vertex*, vector<Element<Dim>*>> vertexElements = BuildVertexElementMap();

		// Sort the vertices by number of elements
		sort(this->Vertices.begin(), this->Vertices.end(), [&vertexElements](Vertex* v1, Vertex* v2) { return vertexElements[v1].size() > vertexElements[v2].size(); });

		for (Vertex* v : this->Vertices)
		{
			vector<Element<Dim>*> elementsAroundV = vertexElements[v];

			if (elementsAroundV.empty()) // the vertex has been discarded (see below)
				continue;

			/*bool vIsOnBoundary = true;
			for (Element<Dim>* e : elementsAroundV)
			{
				if (!e->IsOnBoundary())
				{
					vIsOnBoundary = false;
					break;
				}
			}
			if (vIsOnBoundary)
				continue;*/

			for (Element<Dim>* e : elementsAroundV)
				assert(e->CoarserElement == nullptr);

			coarseMesh->AgglomerateByVertexRemoval(elementsAroundV, v);

			// All the vertices of the agglomerated elements are dicarded for future agglomeration
			for (Element<Dim>* e : elementsAroundV)
			{
				for (Vertex* v2 : e->Vertices())
					vertexElements[v2].clear();
			}
		}
		this->FinalizeCoarsening();
	}


	map<Vertex*, vector<Element<Dim>*>> BuildVertexElementMap()
	{
		// Init the map with an empty list of elements for each vertex
		map<Vertex*, vector<Element<Dim>*>> vertexElements;
		for (Vertex* v : this->Vertices)
			vertexElements.insert({ v, vector<Element<Dim>*>() });

		// Add the elements to the list of their respective vertices
		for (Element<Dim>* e : this->Elements)
		{
			for (Vertex* v : e->Vertices())
				vertexElements[v].push_back(e);
		}

		// TODO: sort the lists in direct order


		return vertexElements;
	}

	/*set<pair<Vertex*, int>> SortByNumberOfElements(map<Vertex*, vector<Element<Dim>*>> vertexElements)
	{
		set<pair<Vertex*, int>> orderedList;

		typename map<Vertex*, vector<Element<Dim>*>>::iterator it;
		for (it = vertexElements.begin(); it != vertexElements.end(); it++)
		{
			Vertex* v = it->first;
			vector<Element<Dim>*> vElements = it->second;

			orderedList.insert(make_pair(v, vElements->size()));
		}

		return orderedList;
	}*/

	void SortInDirectOrder(vector<Element<Dim>*> elementsAroundV, Vertex* v)
	{
		assert(elementsAroundV.size() > 0);
		Element<Dim>* firstElement = elementsAroundV[0];

		// TODO
	}

	//-------------------------------------------------------------//
	// Creates the macro-element obtained after removing a vertex, //
	// adds it to the mesh, and return it.                         //
	//-------------------------------------------------------------//
	Element<Dim>* AgglomerateByVertexRemoval(const vector<Element<Dim>*>& fineElements, Vertex* removedVertex)
	{
		//-----------------------------------------------------------------------------------------------------------//
		// Requirement: fineElements must be sorted in direct (= counter-clockwise) order around the removed vertex. //
		//-----------------------------------------------------------------------------------------------------------//

		// The interior faces of the future macro-element are flagged.
		for (Element<Dim>* e1 : fineElements)
		{
			for (Element<Dim>* e2 : fineElements)
			{
				if (e1 != e2)
				{
					Face<Dim>* interface = e1->CommonFaceWith(e2);
					if (interface != nullptr)
						interface->IsRemovedOnCoarserGrid = true;
				}
			}
		}

		// Separation of the faces we keep and the faces we remove,
		// and the vertices of the macro-element are extracted from those of the faces we keep.
		vector<Vertex*> macroElementVertices;
		set<Face<Dim>*> keptFaces;
		set<Face<Dim>*> removedFaces;
		for (Element<Dim>* e : fineElements)
		{
			for (Face<Dim>* f : e->Faces)
			{
				if (!f->IsRemovedOnCoarserGrid)
				{
					auto ret = keptFaces.insert(f);
					if (ret.second) // if f isn't already in the list keptFaces
					{
						for (Vertex* v : f->Vertices())
						{
							// TODO: Comment s'assurer qu'on les ajoute bien dans l'ordre direct ?
							if (find(macroElementVertices.begin(), macroElementVertices.end(), v) == macroElementVertices.end())
								macroElementVertices.push_back(v);
						}
					}
				}
				else
					removedFaces.insert(f);
			}
		}

		// Creation of the polygonal macro-element
		Element<Dim>* macroElement = CreatePolyhedron(macroElementVertices);
		macroElement->PhysicalPart = fineElements[0]->PhysicalPart;

		// Link between macro-element and fine elements
		for (Element<Dim>* e : fineElements)
		{
			macroElement->FinerElements.push_back(e);
			e->CoarserElement = macroElement;
		}

		// Kept faces are cloned for the coarse mesh and linked to their clones and the macro-element.
		for (Face<Dim>* f : keptFaces)
		{
			CloneAndAddFace(f, macroElement);
		}

		for (Face<Dim>* f : removedFaces)
			macroElement->FinerFacesRemoved.push_back(f);

		// Add macro-element to the list
		this->AddElement(macroElement);

		return macroElement;
	}


	//---------------------------------------------------------------//
	//                                                               //
	//        Coarsening by successive pairwise agglomeration        //
	//                                                               //
	//---------------------------------------------------------------//

	void CoarsenByAgglomerationByPairs(CoarseningStrategy strategy)
	{
		// Intermediary coarsenings by pairs of elements //
		PolyhedralMesh<Dim>* coarseMesh = this;
		for (int i = 0; i < Dim; i++)
		{
			coarseMesh->AggregatePairsOfElements(strategy);
			coarseMesh = dynamic_cast<PolyhedralMesh<Dim>*>(coarseMesh->CoarseMesh);
		}

		RemoveIntermediaryCoarsenings(coarseMesh);
		this->FinalizeCoarsening();
	}

	//----------------------------------------------------------------------------//
	// Relink this' elements and faces to their coarse counterparts in coarseMesh //
	//----------------------------------------------------------------------------//
	void RemoveIntermediaryCoarsenings(PolyhedralMesh<Dim>* coarseMesh)
	{
		for (Element<Dim>* ce : coarseMesh->Elements)
		{
			ce->FinerElements.clear();
			ce->FinerFacesRemoved.clear();
		}

		for (Element<Dim>* fe : this->Elements)
		{
			Element<Dim>* ce = fe->CoarserElement;
			while (ce->CoarserElement)
				ce = ce->CoarserElement;
			fe->CoarserElement = ce;
			ce->FinerElements.push_back(fe);
		}

		for (Face<Dim>* cf : coarseMesh->Faces)
			cf->FinerFaces.clear();

		for (Face<Dim>* ff : this->Faces)
		{
			if (ff->IsRemovedOnCoarserGrid)
			{
				Element<Dim>* ce = ff->Element1->CoarserElement;
				ce->FinerFacesRemoved.push_back(ff);
			}
			else
			{
				Face<Dim>* cf = ff->CoarseFace;
				while (cf->CoarseFace)
					cf = cf->CoarseFace;
				ff->CoarseFace = cf;
				cf->FinerFaces.push_back(ff);
			}
		}

		// Delete the intermediary meshes
		coarseMesh->FineMesh->CoarseMesh = nullptr;
		//delete this->CoarseMesh;

		// Merge faces

		// Link to the coarse mesh
		this->CoarseMesh = coarseMesh;
		coarseMesh->FineMesh = this;
	}

	//-------------------------------------------------------------------------------------------------------------------------//
	// Creates a coarse mesh by aggregates of pairs of elements according to a criterion (closest center or largest interface) //
	//-------------------------------------------------------------------------------------------------------------------------//
	void AggregatePairsOfElements(CoarseningStrategy strategy)
	{
		PolyhedralMesh<Dim>* coarseMesh = new PolyhedralMesh<Dim>();
		this->InitializeCoarsening(coarseMesh);

		for (Element<Dim>* e : this->Elements)
		{
			if (e->CoarserElement)
				continue; // element already aggregated
			
			// Find the neighbour still available for aggregation under a distance/largest interface criterion
			Element<Dim>* neighbourForAggreg = this->FittestAvailableNeighbour(e, strategy);

			if (neighbourForAggreg)
			{
				//-------------------------------------------------------------//
				//           Agglomeration with another fine element           //
				//-------------------------------------------------------------//
				Element<Dim>* macroElement = coarseMesh->AgglomerateFineElements(e, neighbourForAggreg);
				coarseMesh->TryCollapseInterfacesMadeOfMultipleFaces(macroElement);
			}
			else
			{
				// If all the neighbours have already been aggregated, then we aggregate with a macroElement under the same criterion
				Element<Dim>* coarseNeighbourForAggreg = this->FittestCoarseNeighbour(e, strategy);

				if (!coarseNeighbourForAggreg)
					Utils::FatalError("Element cannot be aggregated. Weird...");

				//--------------------------------------------------------//
				//           Agglomeration with a coarse element          //
				//--------------------------------------------------------//
				Element<Dim>* macroElement = coarseMesh->AgglomerateFineElementToCoarse(e, coarseNeighbourForAggreg);
				coarseMesh->TryCollapseInterfacesMadeOfMultipleFaces(macroElement);
			}
		}

		for (Element<Dim>* ce : coarseMesh->Elements)
			coarseMesh->TryCollapseInterfacesMadeOfMultipleFaces(ce);
	}



	//--------------------------------------//
	//                                      //
	//        General useful methods        //
	//                                      //
	//--------------------------------------//

	void RemoveCoarsenedElements(vector<Element<Dim>*>& list)
	{
		vector<Element<Dim>*> uncoarsened;
		uncoarsened.reserve(list.size());
		copy_if(list.begin(), list.end(), back_inserter(uncoarsened), [](Element<Dim>* e) { return !e->CoarserElement; });
		list = uncoarsened;
	}

	vector<Element<Dim>*> AvailableFaceNeighbours(Element<Dim>* elem)
	{
		vector<Element<Dim>*> availableNeighbours;
		for (Element<Dim>* n : elem->NeighboursInSamePhysicalPart())
		{
			if (!n->CoarserElement)
				availableNeighbours.push_back(n);
		}

		return availableNeighbours;
	}

	vector<Element<Dim>*> LockAvailableFaceNeighbours(Element<Dim>* elem, bool doesntFitIfShareReEntrantCorner = false)
	{
		vector<Vertex*> reentrantCorners;
		if (doesntFitIfShareReEntrantCorner)
			reentrantCorners = ReEntrantCorners(elem);

		vector<Element<Dim>*> availableNeighbours;
		for (Element<Dim>* n : elem->NeighboursInSamePhysicalPart())
		{
			bool checkReEntrantCornersOk = !doesntFitIfShareReEntrantCorner || !n->HasAny(reentrantCorners);
			if (!n->CoarserElement && checkReEntrantCornersOk && n->Mutex.try_lock())
			{
				availableNeighbours.push_back(n);
				if (doesntFitIfShareReEntrantCorner)
					reentrantCorners = Utils::Join(reentrantCorners, ReEntrantCorners(n));
			}
		}

		return availableNeighbours;
	}

	vector<Element<Dim>*> AvailableVertexNeighbours(Element<Dim>* elem, map<Vertex*, vector<Element<Dim>*>>& vertexElements)
	{
		set<Element<Dim>*> availableNeighboursSet;
		for (Vertex* v : elem->Vertices())
		{
			auto vertexNeighbours = vertexElements[v];
			for (Element<Dim>* e : vertexNeighbours)
			{
				if (e != elem && !e->CoarserElement)
					availableNeighboursSet.insert(e);
			}
		}

		return vector<Element<Dim>*>(availableNeighboursSet.begin(), availableNeighboursSet.end());
	}

	Element<Dim>* FittestAvailableNeighbour(Element<Dim>* e, CoarseningStrategy strategy)
	{
		// Find the neighbour still available for aggregation under a distance/largest interface criterion
		Element<Dim>* neighbourForAggreg = nullptr;
		double closestDistance = -1;
		double largestInterface = -1;

		for (Face<Dim>* f : e->Faces)
		{
			if (f->IsDomainBoundary || f->IsRemovedOnCoarserGrid)
				continue;

			Element<Dim>* neighbour = f->GetNeighbour(e);
			if (neighbour->CoarserElement) // already aggregated
				continue;

			// Init choice criterion
			double distance = 0;
			double interfaceMeasure = 0;
			if (strategy == CoarseningStrategy::AgglomerationCoarseningByClosestCenter)
				distance = Vect<Dim>(e->Center(), neighbour->Center()).norm();
			else if (strategy == CoarseningStrategy::AgglomerationCoarseningByClosestFace)
				distance = Vect<Dim>(e->Center(), f->Center()).norm();
			else if (strategy == CoarseningStrategy::AgglomerationCoarseningByLargestInterface)
			{
				for (Face<Dim>* fInterface : e->Faces)
				{
					if (fInterface->IsDomainBoundary || fInterface->IsRemovedOnCoarserGrid || fInterface->GetNeighbour(e) != neighbour)
						continue;
					interfaceMeasure += fInterface->Measure();
				}
			}

			if (!neighbourForAggreg)
			{
				neighbourForAggreg = neighbour;
				closestDistance = distance;
				largestInterface = interfaceMeasure;
			}
			else
			{
				if (strategy == CoarseningStrategy::AgglomerationCoarseningByClosestCenter && distance < closestDistance)
				{
					neighbourForAggreg = neighbour;
					closestDistance = distance;
				}
				else if (strategy == CoarseningStrategy::AgglomerationCoarseningByClosestFace && distance < closestDistance)
				{
					neighbourForAggreg = neighbour;
					closestDistance = distance;
				}
				else if (strategy == CoarseningStrategy::AgglomerationCoarseningByLargestInterface && interfaceMeasure > largestInterface)
				{
					neighbourForAggreg = neighbour;
					largestInterface = interfaceMeasure;
				}
			}
		}

		return neighbourForAggreg;
	}

	Element<Dim>* FittestCoarseNeighbour(Element<Dim>* e, CoarseningStrategy strategy)
	{
		bool mustCheckReEntrantCorners = Utils::ProgramArgs.Solver.MG.ReEntrantCornerManagement != ReEntrantCornerMgmt::Disabled;
		vector<Vertex*> reentrantCorners = ReEntrantCorners(e);

		Element<Dim>* coarseNeighbourForAggreg = nullptr;
		double closestDistance = -1;
		double largestInterface = -1;
		for (Face<Dim>* f : e->Faces)
		{
			if (f->IsDomainBoundary)
				continue;

			Element<Dim>* macroNeighbour = f->GetNeighbour(e)->CoarserElement;
			if (!macroNeighbour || !macroNeighbour->IsInSamePhysicalPartAs(e))
				continue;

			if (mustCheckReEntrantCorners && macroNeighbour->HasAny(reentrantCorners))
				continue;

			double distance = 0;
			double interfaceMeasure = 0;
			if (strategy == CoarseningStrategy::AgglomerationCoarseningByClosestCenter)
				distance = Vect<Dim>(e->Center(), macroNeighbour->Center()).norm();
			else if (strategy == CoarseningStrategy::AgglomerationCoarseningByClosestFace)
				distance = Vect<Dim>(e->Center(), f->Center()).norm();
			else if (strategy == CoarseningStrategy::AgglomerationCoarseningByLargestInterface)
			{
				for (Face<Dim>* fInterface : e->Faces)
				{
					if (fInterface->IsDomainBoundary || fInterface->IsRemovedOnCoarserGrid || fInterface->GetNeighbour(e)->CoarserElement != macroNeighbour)
						continue;
					interfaceMeasure += fInterface->Measure();
				}
			}

			if (!coarseNeighbourForAggreg)
			{
				coarseNeighbourForAggreg = macroNeighbour;
				closestDistance = distance;
				largestInterface = interfaceMeasure;
			}
			else
			{
				if (strategy == CoarseningStrategy::AgglomerationCoarseningByClosestCenter && distance < closestDistance)
				{
					coarseNeighbourForAggreg = macroNeighbour;
					closestDistance = distance;
				}
				else if (strategy == CoarseningStrategy::AgglomerationCoarseningByClosestFace && distance < closestDistance)
				{
					coarseNeighbourForAggreg = macroNeighbour;
					closestDistance = distance;
				}
				else if (strategy == CoarseningStrategy::AgglomerationCoarseningByLargestInterface && interfaceMeasure > largestInterface)
				{
					coarseNeighbourForAggreg = macroNeighbour;
					largestInterface = interfaceMeasure;
				}
			}
		}

		return coarseNeighbourForAggreg;
	}

	vector<Vertex*> ReEntrantCorners(Element<Dim>* e)
	{
		vector<Vertex*> reentrantCorners;
		if (Utils::ProgramArgs.Solver.MG.ReEntrantCornerManagement != ReEntrantCornerMgmt::Disabled)
		{
			auto it = this->_reEntrantCorners.find(e->PhysicalPart);
			if (it != this->_reEntrantCorners.end())
			{
				for (Vertex* corner : it->second)
				{
					if (e->HasVertex(corner))
						reentrantCorners.push_back(corner);
				}
			}
		}
		return reentrantCorners;
	}


	//--------------------------------------------------------------------//
	// Creates the macro-element obtained by agglomeration of 2 elements, //
	// adds it to the mesh, and return it.                                //
	//--------------------------------------------------------------------//
	Element<Dim>* AgglomerateFineElements(Element<Dim>* e1, Element<Dim>* e2)
	{
		assert(e1->IsInSamePhysicalPartAs(e2));

		// Collect the faces interfacing these two elements //
		vector<Face<Dim>*> facesToRemove;
		vector<Face<Dim>*> facesToClone;
		for (Face<Dim>* f : e1->Faces)
		{
			if (f->GetNeighbour(e1) == e2)
				facesToRemove.push_back(f);
			else
				facesToClone.push_back(f);
		}
		for (Face<Dim>* f : e2->Faces)
		{
			if (f->GetNeighbour(e2) != e1)
				facesToClone.push_back(f);
		}

		// Creation of the polygonal macro-element
		Element<Dim>* coarseElement = CreateMacroElement(e1, e2, facesToRemove);
		coarseElement->PhysicalPart = e1->PhysicalPart;

		// Links between the coarse element and the fine elements
		coarseElement->FinerElements.push_back(e1);
		coarseElement->FinerElements.push_back(e2);
		e1->CoarserElement = coarseElement;
		e2->CoarserElement = coarseElement;

		for (Face<Dim>* f : facesToRemove)
		{
			f->IsRemovedOnCoarserGrid = true;
			coarseElement->FinerFacesRemoved.push_back(f);
		}

		// Add coarse element to the list
		this->AddElement(coarseElement);

		// Kept faces are cloned for the coarse mesh and linked to their clones and the coarse element.
		for (Face<Dim>* f : facesToClone)
			this->CloneAndAddFace(f, coarseElement);

		assert(coarseElement->Faces.size() > 2);
		if (Dim == 2)
			assert(coarseElement->Faces.size() == coarseElement->Vertices().size());

		return coarseElement;
	}

	Element<Dim>* AgglomerateFineElements(const vector<Element<Dim>*>& fineElements, bool checkNoSharedReEntrantCorner = false)
	{
		//assert(fineElements.size() > 1 || fineElements[0]->IsOnBoundary());
		for (Element<Dim>* e : fineElements)
		{
			assert(!e->CoarserElement);
			assert(e->IsInSamePhysicalPartAs(fineElements[0]) && "Agglomeration of elements from different physical groups is not allowed!");
		}

		// Check re-entrant corners
		if (checkNoSharedReEntrantCorner)
		{
			for (int i = 0; i < fineElements.size(); i++)
			{
				vector<Vertex*> reentrantCorners = ReEntrantCorners(fineElements[i]);
				for (int j = i + 1; j < fineElements.size(); j++)
				{
					if (fineElements[j]->HasAny(reentrantCorners))
					{
						MatlabScript s;
						s.OpenFigure();
						for (Element<Dim>* e : fineElements)
							s.PlotPolygonEdges(e->Shape()->Vertices(), "k");
						fineElements[i]->ExportToMatlab("r", true);
						fineElements[j]->ExportToMatlab("r", true);
						for (Vertex* rc : reentrantCorners)
							s.PlotPoint(*rc, "b+");
						Utils::FatalError("Agglomeration of elements sharing a reentrant corner is not allowed.");
					}
				}
			}
		}

		Agglo<Dim> agglo(fineElements);
		if (!agglo.Success)
			return nullptr;

		Element<Dim>* coarseElement = CreatePolyhedron(agglo.Vertices());
		coarseElement->PhysicalPart = fineElements[0]->PhysicalPart;

		assert(coarseElement->Vertices().size() > 2);

		for (Element<Dim>* e : fineElements)
		{
			e->CoarserElement = coarseElement;
			coarseElement->FinerElements.push_back(e);
			if (e->HasReEntrantCorner)
				coarseElement->HasReEntrantCorner = true;
		}

		for (Face<Dim>* f : agglo.RemovedFaces())
		{
			f->IsRemovedOnCoarserGrid = true;
			coarseElement->FinerFacesRemoved.push_back(f);
		}

		// Add coarse element to the list
		this->AddElement(coarseElement);

		// Kept faces are cloned for the coarse mesh and linked to their clones and the coarse element.
		for (Face<Dim>* f : agglo.Faces)
			this->CloneAndAddFace(f, coarseElement);

		assert(coarseElement->Faces.size() > 2);
		if (Dim == 2)
			assert(coarseElement->Faces.size() == coarseElement->Vertices().size());

		return coarseElement;
	}

	Element<Dim>* AgglomerateFineElementToCoarse(Element<Dim>* fineElement, Element<Dim>* coarseElement, bool checkNoSharedReEntrantCorner = false)
	{
		assert(fineElement->IsInSamePhysicalPartAs(coarseElement) && "Agglomeration of elements from different physical groups is not allowed!");
		
		// Collect the fine faces interfacing these two elements //
		vector<Face<Dim>*> facesToRemove;
		vector<Face<Dim>*> facesToClone;
		vector<Face<Dim>*> coarseFacesToKeep;
		for (Face<Dim>* ff : fineElement->Faces)
		{
			if (!ff->IsDomainBoundary && ff->CoarseFace && coarseElement->HasFace(ff->CoarseFace))
				facesToRemove.push_back(ff);
			else
				facesToClone.push_back(ff);
		}
		for (Face<Dim>* cf : coarseElement->Faces)
		{
			bool mustBeKept = true;
			for (Face<Dim>* ff : cf->FinerFaces)
			{
				for (Face<Dim>* fToRemove : facesToRemove)
				{
					if (ff == fToRemove)
					{
						mustBeKept = false;
						break;
					}
				}
				if (!mustBeKept)
					break;
			}
			if (mustBeKept)
				coarseFacesToKeep.push_back(cf);
		}

		// Check re-entrant corners
		if (checkNoSharedReEntrantCorner)
		{
			vector<Vertex*> reentrantCorners = ReEntrantCorners(fineElement);
			for (Face<Dim>* f : facesToRemove)
			{
				if (f->HasAny(reentrantCorners))
					Utils::FatalError("Agglomeration of elements sharing a reentrant corner is not allowed.");
			}
		}
		
		// Creation of the polygonal macro-element
		Element<Dim>* newCoarseElement;
		try
		{
			newCoarseElement = CreateMacroElement(fineElement, coarseElement, facesToRemove);
		}
		catch (exception* e)
		{
			string error(e->what());
			Utils::Error("An exception occured during the agglomeration of a fine and coarse elements: " + error);
			cout << "% ------------------- Fine element:" << endl;
			fineElement->ExportToMatlab("r");
			cout << "% ------------------- Coarse element:" << endl;
			coarseElement->ExportToMatlab("b");
			throw e;
		}
		newCoarseElement->Mutex.lock();
		newCoarseElement->Id = this->NewElementId();
		newCoarseElement->PhysicalPart = coarseElement->PhysicalPart;
		newCoarseElement->HasReEntrantCorner = coarseElement->HasReEntrantCorner;

		// Links between the macro element and the fine elements
		for (Element<Dim>* fe : coarseElement->FinerElements)
		{
			newCoarseElement->FinerElements.push_back(fe);
			fe->CoarserElement = newCoarseElement;
		}
		newCoarseElement->FinerElements.push_back(fineElement);
		fineElement->CoarserElement = newCoarseElement;
		if (fineElement->HasReEntrantCorner)
			newCoarseElement->HasReEntrantCorner = true;

		for (Face<Dim>* ff : coarseElement->FinerFacesRemoved)
			newCoarseElement->FinerFacesRemoved.push_back(ff);

		for (Face<Dim>* cf : coarseFacesToKeep)
		{
			if (cf->Element1 == coarseElement)
				cf->Element1 = newCoarseElement;
			else
				cf->Element2 = newCoarseElement;
			newCoarseElement->Faces.push_back(cf);
		}

		for (Face<Dim>* ff : facesToRemove)
		{
			ff->IsRemovedOnCoarserGrid = true;
			newCoarseElement->FinerFacesRemoved.push_back(ff);
			if (ff->CoarseFace->FinerFaces.size() > 1 && !ff->CoarseFace->IsDeleted)
				RemoveFromCoarseFacesNotToCross(ff->CoarseFace);
			this->RemoveFace(ff->CoarseFace, false);
			ff->CoarseFace = nullptr;
		}

		// Delete the old coarse element and replace it with the new one
		this->ReplaceAndDeleteElement(coarseElement, newCoarseElement);

		// Kept faces are cloned for the coarse mesh and linked to their clones and the coarse element.
		for (Face<Dim>* ff : facesToClone)
			this->CloneAndAddFace(ff, newCoarseElement);

		assert(newCoarseElement->Faces.size() > 2);
		if (Dim == 2)
			assert(newCoarseElement->Faces.size() == newCoarseElement->Vertices().size());

		newCoarseElement->Mutex.unlock();
		return newCoarseElement;
	}

	Element<Dim>* Agglomerate(Element<Dim>* e1, Element<Dim>* e2)
	{
		vector<Element<Dim>*> v{ e1, e2 };
		return Agglomerate(v);
	}

	Element<Dim>* Agglomerate(const vector<Element<Dim>*>& elements)
	{
		Agglo<Dim> agglo(elements);
		Element<Dim>* newElement = CreatePolyhedron(agglo.Vertices());
		newElement->Id = this->NewElementId();
		newElement->PhysicalPart = elements[0]->PhysicalPart;

		for (Element<Dim>* e : elements)
		{
			assert(e->IsInSamePhysicalPartAs(elements[0]) && "Agglomeration of elements from different physical groups is not allowed!");
			for (Face<Dim>* f : e->Faces)
			{
				if (f->Element1 == e)
					f->Element1 = newElement;
				if (f->Element2 == e)
					f->Element2 = newElement;
			}
			for (Element<Dim>* fe : e->FinerElements)
			{
				fe->CoarserElement = newElement;
				newElement->FinerElements.push_back(fe);
			}
			for (Face<Dim>* ff : e->FinerFacesRemoved)
				newElement->FinerFacesRemoved.push_back(ff);
		}

		for (Face<Dim>* f : agglo.RemovedFaces())
		{
			for (Face<Dim>* ff : f->FinerFaces)
			{
				ff->IsRemovedOnCoarserGrid = true;
				newElement->FinerFacesRemoved.push_back(ff);
			}
			// If it's a collapsed face, we also need to remove its references in CoarseFacesNotToCross
			if (f->FinerFaces.size() > 1)
				RemoveFromCoarseFacesNotToCross(f);
			this->RemoveFace(f, false);
		}

		newElement->Faces = agglo.Faces;

		// Replace the first element with the new one
		this->ReplaceAndDeleteElement(elements[0], newElement);
		// ... and remove the others
		for (int i = 1; i < elements.size(); i++)
			this->RemoveElement(elements[i]);

		assert(newElement->Vertices().size() > 2);
		assert(newElement->Faces.size() > 2);
		if (Dim == 2)
			assert(newElement->Faces.size() == newElement->Vertices().size());

		return newElement;
	}

	FaceCollapsingStatus AgglomerateFineFaces(Face<Dim>* f1, Face<Dim>* f2)
	{
		// f1 and f2 are one level higher than 'this'
		return TryCollapse({ f1->CoarseFace, f2->CoarseFace });
	}

	FaceCollapsingStatus TryCollapseByPairs(const vector<Face<Dim>*>& interfaceFaces)
	{
		set<Face<Dim>*> faces;
		for (Face<Dim>* f : interfaceFaces)
		{
			assert(f->FinerFaces.size() < 3);
			if (f->FinerFaces.size() == 1)
				faces.insert(f);
		}

		FaceCollapsingStatus returnStatus = FaceCollapsingStatus::NotEnoughFaces;
		while (faces.size() > 1)
		{
			auto it1 = faces.begin();
			Face<Dim>* f1 = *it1;
			auto it2 = it1;
			it2++;
			while (it2 != faces.end())
			{
				Face<Dim>* f2 = *it2;
				FaceCollapsingStatus status = TryCollapse({ f1, f2 });
				if (status == FaceCollapsingStatus::Ok)
				{
					returnStatus = FaceCollapsingStatus::Ok;
					faces.erase(it2);
					break;
				}
				else if (returnStatus != FaceCollapsingStatus::Ok && status != FaceCollapsingStatus::InterfaceHasHoles)
					returnStatus = status;
				it2++;
			}
			faces.erase(it1);
		}
		return returnStatus;
	}

	FaceCollapsingStatus TryCollapse(const vector<Face<Dim>*>& faces)
	{
		assert(faces.size() > 1);
		// 'faces' must be at the same level as 'this'

		Interface<Dim> interf(faces);

		assert(!this->_geometricVertices.empty());
		for (Vertex* v : interf.InteriorVertices())
		{
			if (this->_geometricVertices.find(v) != this->_geometricVertices.end())
				return FaceCollapsingStatus::GeometricErosion;
		}

		FaceCollapsingStatus status = interf.AnalyzeCollapsing();
		if (status != FaceCollapsingStatus::Ok)
			return status;

		Collapse(interf);
		return FaceCollapsingStatus::Ok;
	}

	void Collapse(Interface<Dim>& interf)
	{
		// The faces and mergedFace must be at the same level.

		vector<Face<Dim>*> faces = interf.Faces();
		Face<Dim>* mergedFace = interf.CollapsedFace();

		mergedFace->IsDomainBoundary = faces[0]->IsDomainBoundary;
		for (Face<Dim>* f : faces)
		{
			if (faces[0]->BoundaryPart)
				assert(f->BoundaryPart->Id == faces[0]->BoundaryPart->Id && "To be agglomerated, all faces must have the same BoundaryPart.");
			else
				assert(!f->BoundaryPart && "To be agglomerated, all faces must have the same BoundaryPart.");
		}
		mergedFace->BoundaryPart = faces[0]->BoundaryPart;

		Element<Dim>* elem1 = faces[0]->Element1;
		Element<Dim>* elem2 = faces[0]->Element2;
		mergedFace->Element1 = elem1;
		mergedFace->Element2 = elem2;

		// Intergrid links
		for (Face<Dim>* f : faces)
		{
			for (Face<Dim>* ff : f->FinerFaces)
			{
				mergedFace->FinerFaces.push_back(ff);
				ff->CoarseFace = mergedFace;
			}
		}

		// Replace the faces in the elements (reshape the polygon, update the vertices, transfer finer elements to the neighbours if needed)
		Collapse(elem1, interf);
		if (elem2)
			Collapse(elem2, interf);

		// Remove the faces from the mesh
		for (Face<Dim>* f : faces)
		{
			// If it's a collapsed face, we also need to remove its references in CoarseFacesNotToCross
			if (f->FinerFaces.size() > 1)
				RemoveFromCoarseFacesNotToCross(f);
			this->RemoveFace(f, false);
		}

		// Add mergedFace to the list
		this->AddFace(mergedFace, false);
	}

	static void RemoveFromCoarseFacesNotToCross(Face<Dim>* cf)
	{
		vector<Element<Dim>*> fineElements = cf->Element1->FinerElements;
		if (cf->Element2)
			fineElements = Utils::Join(fineElements, cf->Element2->FinerElements);
		for (Element<Dim>* fe : fineElements)
		{
			auto it = fe->CoarseFacesNotToCross.find(cf);
			if (it != fe->CoarseFacesNotToCross.end())
				fe->CoarseFacesNotToCross.erase(it);
		}
	}

	// Replace faces with their agglomeration //
	void Collapse(Element<Dim>* e, Interface<Dim> interf)
	{
		vector<Face<Dim>*> faces = interf.Faces();
		Face<Dim>* collapsedFace = interf.CollapsedFace();

		vector<Face<Dim>*> currentFaces = e->Faces;
		e->Faces.clear();
		for (Face<Dim>* f : currentFaces)
		{
			if (!f->IsIn(faces))
				e->Faces.push_back(f);
		}
		e->Faces.push_back(collapsedFace);

		e->RemoveIntersections(interf.InteriorVertices());

		// Some fine elements associated to this coarse one might now be more overlapping the other coarse elements.
		// So we change the association.
		if (!collapsedFace->IsDomainBoundary)
		{
			Element<Dim>* neighbour = collapsedFace->GetNeighbour(e);

			auto it = e->FinerElements.begin();
			while (it != e->FinerElements.end())
			{
				Element<Dim>* fe = *it;
				if (fe->IsInSamePhysicalPartAs(neighbour) && (fe->HasAny(collapsedFace->FinerFaces) || fe->HasAny(interf.InteriorVertices())))
				{
					fe->IsFullyEmbeddedInCoarseElement = false;
					fe->CoarseFacesNotToCross.insert(collapsedFace);
					DomPoint p = fe->InteriorPoint();
					if (!e->Contains(p) && !collapsedFace->Contains(p) && neighbour->Contains(p))
					{
						/*cout << "%---------- this " << endl;
						this->ExportToMatlab("b");
						cout << "% neighbour " << endl;
						neighbour->ExportToMatlab("m");
						cout << "% transfered finer element (b --> m) " << endl;
						fe->ExportToMatlab("r");
						MatlabScript s;
						s.PlotPoint(p, "rx");*/

						// Transfer to the neighbour
						fe->CoarserElement = neighbour;
						neighbour->FinerElements.push_back(fe);
						it = e->FinerElements.erase(it);
						continue;
					}
				}
				it++;
			}
		}

	}

	void CloneAndAddFace(Face<Dim>* f, Element<Dim>* macroElement)
	{
		f->Mutex.lock();
		if (!f->CoarseFace) // this face has not already been created by a previous aggregation
		{
			// Clone face
			BigNumber faceNumber = this->Faces.size();
			Face<Dim>* copy = f->CreateSameGeometricFace(faceNumber, macroElement);

			// Associate
			f->CoarseFace = copy;
			copy->FinerFaces.push_back(f);

			// Add copied face to the list of coarse faces
			this->AddFace(copy, false);
		}
		else
		{
			if (!f->CoarseFace->Element2)
				f->CoarseFace->Element2 = macroElement;

			assert(f->CoarseFace->Element1 != f->CoarseFace->Element2);
		}
		f->Mutex.unlock();
		macroElement->AddFace(f->CoarseFace);
	}

public:
	void SetOverlappingFineElements() override
	{
		CoarseningStrategy stgy = this->ComesFrom.CS;
		if (stgy == CoarseningStrategy::None)
			stgy = this->FineMesh->ComesFrom.CS;

		ElementParallelLoop<Dim> parallelLoop(this->Elements);
		parallelLoop.Execute([stgy](Element<Dim>* ce)
			{
				SetOverlappingFineElements(ce, stgy);
				ce->InitOverlappingElementsLocalNumbering();
			});
	}

private:
	static void SetOverlappingFineElements(Element<Dim>* ce, CoarseningStrategy stgy)
	{
		set<Element<Dim>*> tested;
		double overlappingMeasure = 0;
		bool stop = false;
		for (auto fe : ce->FinerElements)
		{
			assert(fe->PhysicalPart == ce->PhysicalPart);
			CheckOverlapping(ce, fe, tested, overlappingMeasure, stop);
		}

		if (stop)
			return;

		if (!Utils::BuildsNestedMeshHierarchy(stgy))
		{
			for (Element<Dim>* neighbour : ce->VertexNeighbours())
			{
				for (auto fe1 : neighbour->FinerElements)
				{
					CheckOverlapping(ce, fe1, tested, overlappingMeasure, stop);
					if (stop)
						return;
					for (Element<Dim>* fe2 : fe1->VertexNeighbours())
					{
						CheckOverlapping(ce, fe2, tested, overlappingMeasure, stop);
						if (stop)
							return;
						for (Element<Dim>* fe3 : fe2->VertexNeighbours())
						{
							CheckOverlapping(ce, fe3, tested, overlappingMeasure, stop);
							if (stop)
								return;
						}
					}
				}
			}
		}

		if (!ce->IsAtPhysicalPartBoundary())
			Utils::FatalError("error in the computation of the overlapping elements: a coarse element is not fully overlapped by finer ones.");
	}

	static void CheckOverlapping(Element<Dim>* ce, Element<Dim>* fe, set<Element<Dim>*>& tested, double& overlappingMeasure, bool& stop)
	{
		if (ce->IsInSamePhysicalPartAs(fe) && tested.find(fe) == tested.end())
		{
			vector<PhysicalShape<Dim>*> intersection = Intersection(ce, fe);
			if (!intersection.empty())
			{
				ce->OverlappingFineElements.insert({ fe , intersection });
				for (PhysicalShape<Dim>* intersec : intersection)
					overlappingMeasure += intersec->Measure();
				if (abs(overlappingMeasure - ce->Measure()) < Utils::Eps * ce->Measure())
					stop = true;
			}
		}
		tested.insert(fe);
	}

public:
	void SetOverlappingFineElementsSubTriangles() override
	{
		CoarseningStrategy stgy = this->ComesFrom.CS;
		if (stgy == CoarseningStrategy::None)
			stgy = this->FineMesh->ComesFrom.CS;

		ElementParallelLoop<Dim> parallelLoopFine(this->FineMesh->Elements);
		parallelLoopFine.Execute([stgy](Element<Dim>* fe)
			{
				if (!Utils::BuildsNestedMeshHierarchy(stgy) && !fe->IsFullyEmbeddedInCoarseElement)
				{
					if (stgy == CoarseningStrategy::IndependentRemeshing)
						fe->Refine();
					else
						fe->RefineWithoutCoarseOverlap();

					/*if (fe->CoarserElement->Number == 16 && fe->Number == 164) // (16,164) (38, 152) (38, 40)
					{
						MatlabScript s;
						s.Comment("--------------------COARSE-------------------");
						s.PlotPolygonEdges(fe->CoarserElement->Shape()->Vertices(), "k", 3);
						s.Comment("FINE " + to_string(fe->Number));
						for (auto ss : fe->Shape()->RefinedShapes())
							s.PlotPolygon(ss->Vertices(), "r", "-");

						fe->RefineWithoutCoarseOverlap();
					}*/
				}

				SetOverlappingFineElementsSubTriangles(fe, stgy);
			});

		ElementParallelLoop<Dim> parallelLoop(this->Elements);
		parallelLoop.Execute([](Element<Dim>* ce)
			{
				ce->InitOverlappingElementsLocalNumbering();
			});
	}

private:
	static void SetOverlappingFineElementsSubTriangles(Element<Dim>* fe, CoarseningStrategy stgy)
	{
		if (Utils::BuildsNestedMeshHierarchy(stgy) || fe->IsFullyEmbeddedInCoarseElement)
		{
			Element<Dim>* ce = fe->CoarserElement;
			ce->Mutex.lock();
			ce->OverlappingFineElements.insert({ fe , {fe->Shape()} });
			ce->Mutex.unlock();
			return;
		}

		for (PhysicalShape<Dim>* subShape : fe->Shape()->RefinedShapes())
		{
			assert(fe->IsInSamePhysicalPartAs(fe->CoarserElement));
			vector<Element<Dim>*> coarseCandidates = fe->CoarserElement->ThisAndVertexNeighbours();
			bool found = false;
			for (Element<Dim>* ce : coarseCandidates)
			{
				if (ce->Contains(subShape->Center()))
				{
					ce->Mutex.lock();
					auto it = ce->OverlappingFineElements.find(fe);
					if (it == ce->OverlappingFineElements.end())
						ce->OverlappingFineElements.insert({ fe , {subShape} });
					else
						it->second.push_back(subShape);
					ce->Mutex.unlock();
					found = true;
					break;
				}
			}
			if (!found)
			{
#ifndef NDEBUG
				MatlabScript s;
				s.Comment("--------- coarse element");
				fe->CoarserElement->ExportToMatlab("y", true);
				for (Element<Dim>* ce : fe->CoarserElement->ThisAndVertexNeighbours(true))
				{
					s.Comment("--------- coarse candidate");
					ce->ExportToMatlab("b", true);
				}
				s.Comment("--------- fine element " + to_string(fe->Number));
				fe->ExportToMatlab("g");
				s.Comment("--------- subshape");
				subShape->ExportToMatlab("r");
#endif // !NDEBUG
				if (!fe->IsAtPhysicalPartBoundary() && !fe->IsOnBoundary())
					Utils::Warning("A fine subshape's center that is not a domain/physical region boundary is not in any coarse element");
				else if (fe->IsAtPhysicalPartBoundary())
					Utils::Warning("A fine subshape's center (at physical boundary) is not in any coarse element of same physical part");

				// Affectation to the closest coarse element
				Element<Dim>* closestCoarse = nullptr;
				double minDistance = -1;
				for (Element<Dim>* ce : coarseCandidates)
				{
					double distance = Vect<Dim>(subShape->Center(), ce->Center()).norm();
					if (!closestCoarse || distance < minDistance)
					{
						closestCoarse = ce;
						minDistance = distance;
					}
				}
				closestCoarse->Mutex.lock();
				auto it = closestCoarse->OverlappingFineElements.find(fe);
				if (it == closestCoarse->OverlappingFineElements.end())
					closestCoarse->OverlappingFineElements.insert({ fe , {subShape} });
				else
					it->second.push_back(subShape);
				closestCoarse->Mutex.unlock();
			}
		}
	}

private:
	static CGAL::Polygon_2<exactKernel> ConvertToCGALPolygon(Element<Dim>* e)
	{
		assert(Dim == 2);
		CGAL::Polygon_2<exactKernel> cgalPoly;
		const Polygon* poly = dynamic_cast<const Polygon*>(e->Shape());
		if (poly)
			cgalPoly = poly->CGALPolygon();
		else
			cgalPoly = CGALWrapper::CreatePolygon<exactKernel>(e->Shape()->Vertices());
		return cgalPoly;
	}

public:
	virtual void SanityCheck() override
	{
		Mesh<Dim>::SanityCheck();

		/*if (this->ComesFrom.CS == CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours)
		{
			for (Element<Dim>* e : this->Elements)
			{
				for (Element<Dim>* n : e->NeighboursInSamePhysicalPart())
				{
					auto faces = e->InterfaceWith(n);
					if (faces.size() != 1)
					{
						FaceCollapsingStatus status = TryCollapseInterfaceBetween(e, n);
						if (status == FaceCollapsingStatus::Ok)
						{
							e->ExportToMatlab("r");
							n->ExportToMatlab("m");
							assert(faces.size() == 1 && "After face collapsing, neighbours in the same physical part should have only one face in common.");
						}
					}
				}
			}
		}*/
	}

	virtual ~PolyhedralMesh()
	{
		if (_quadrilateralElements.empty() && _triangularElements.empty() && _tetrahedralElements.empty() && _parallelepipedElements.empty())
		{
			for (size_t i = 0; i < this->Elements.size(); ++i)
				delete this->Elements[i];
			this->Elements.clear();
		}

		if (_edgeFaces.empty() && _triangularFaces.empty())
		{
			for (size_t i = 0; i < this->Faces.size(); ++i)
				delete this->Faces[i];
			this->Faces.clear();
		}
	}

	//-----------------------------------------------//
	// Dim-specific functions (implementation below) //
	//-----------------------------------------------//
	Element<Dim>* CreatePolyhedron(vector<Vertex*> vertices) { return nullptr; }
	Element<Dim>* CreateMacroElement(Element<Dim>* e1, Element<Dim>* e2, const vector<Face<Dim>*>& facesToRemove) { return nullptr; }
	Face<Dim>* CreateMacroFace(Face<Dim>* f1, Face<Dim>* f2, Vertex* vertexToRemove) { return nullptr; }
	static vector<PhysicalShape<Dim>*> Intersection(Element<Dim>* e1, Element<Dim>* e2) { assert(false); }
	void FaceCoarsening() { assert(false); };
};


template<>
Element<2>* PolyhedralMesh<2>::CreatePolyhedron(vector<Vertex*> vertices)
{
	BigNumber elementNumber = this->Elements.size();
	PolygonalElement* macroElement = new PolygonalElement(elementNumber, vertices, false);
	return macroElement;
}

template<>
Element<3>* PolyhedralMesh<3>::CreatePolyhedron(vector<Vertex*> vertices)
{
	assert(false && "Not yet implemented");
}

template<>
Element<2>* PolyhedralMesh<2>::CreateMacroElement(Element<2>* e1, Element<2>* e2, const vector<Face<2>*>& facesToRemove)
{
	PolygonalElement* macroElement = new PolygonalElement(0, e1, e2, facesToRemove, false);
	return macroElement;
}

template<>
Element<3>* PolyhedralMesh<3>::CreateMacroElement(Element<3>* e1, Element<3>* e2, const vector<Face<3>*>& facesToRemove)
{
	Utils::FatalError("Not implemented in 3D.");
	return nullptr;
}

template<>
Face<2>* PolyhedralMesh<2>::CreateMacroFace(Face<2>* f1, Face<2>* f2, Vertex* vertexToRemove)
{
	Vertex* v1 = f1->Vertices()[0] == vertexToRemove ? f1->Vertices()[1] : f1->Vertices()[0];
	Vertex* v2 = f2->Vertices()[0] == vertexToRemove ? f2->Vertices()[1] : f2->Vertices()[0];
	Edge* macroFace = new Edge(0, v1, v2);
	return macroFace;
}

template<>
Face<3>* PolyhedralMesh<3>::CreateMacroFace(Face<3>* f1, Face<3>* f2, Vertex* vertexToRemove)
{
	Utils::FatalError("Not implemented in 3D.");
	return nullptr;
}

template<>
vector<PhysicalShape<2>*> PolyhedralMesh<2>::Intersection(Element<2>* e1, Element<2>* e2)
{
	CGAL::Polygon_2<exactKernel> cgalPoly1 = ConvertToCGALPolygon(e1);
	CGAL::Polygon_2<exactKernel> cgalPoly2 = ConvertToCGALPolygon(e2);

	vector<CGAL::Polygon_2<exactKernel>> inters = CGALWrapper::Intersection(cgalPoly1, cgalPoly2);

	vector<PhysicalShape<2>*> intersectionPolygons;
	for (const CGAL::Polygon_2<exactKernel>& p : inters)
		intersectionPolygons.push_back(new Polygon(p, true));
	return intersectionPolygons;
}

template <>
void PolyhedralMesh<2>::FaceCoarsening()
{
	PolyhedralMesh<2>* coarseSkeleton = new PolyhedralMesh<2>();
	this->InitializeCoarsening(coarseSkeleton);
	coarseSkeleton->ComesFrom.CS = CoarseningStrategy::FaceCoarsening;

	// Copy all vertices
	map<size_t, MeshVertex<2>*> verticesByNumber;
	for (Vertex* v : this->Vertices)
	{
		MeshVertex<2>* coarseV = new MeshVertex<2>(*v);
		coarseSkeleton->Vertices.push_back(coarseV);
		verticesByNumber.insert({ coarseV->Number, coarseV });
	}

	BigNumber faceNumber = 0;

	for (Face<2>* face1 : this->Faces)
	{
		if (face1->CoarseFace)
			continue;

		Edge* edge1 = dynamic_cast<Edge*>(face1);
		CartesianEdge* cedge1 = dynamic_cast<CartesianEdge*>(face1);
		for (Vertex* v : face1->Vertices())
		{
			MeshVertex<2>* vertex = static_cast<MeshVertex<2>*>(v);
			if (!vertex)
				Utils::FatalError("Face coarsening not implemented on this mesh");

			for (Face<2>* face2 : vertex->Faces)
			{
				if (face1 == face2 || face2->CoarseFace)
					continue;

				Edge* edge2 = dynamic_cast<Edge*>(face2);
				CartesianEdge* cedge2 = dynamic_cast<CartesianEdge*>(face2);

				// Compute colinearity
				Vertex* otherVertexInEdge1;
				Vertex* otherVertexInEdge2;
				if (edge1)
				{
					otherVertexInEdge1 = edge1->Vertex1() == v ? edge1->Vertex2() : edge1->Vertex1();
					otherVertexInEdge2 = edge2->Vertex1() == v ? edge2->Vertex2() : edge2->Vertex1();
				}
				else
				{
					vector<Vertex*> verticesEdge1 = face1->Vertices();
					vector<Vertex*> verticesEdge2 = face2->Vertices();
					otherVertexInEdge1 = verticesEdge1[0] == v ? verticesEdge1[1] : verticesEdge1[0];
					otherVertexInEdge2 = verticesEdge2[0] == v ? verticesEdge2[1] : verticesEdge2[0];
				}

				MeshVertex<2>* coarseV = verticesByNumber.at(v->Number);
				MeshVertex<2>* coarseOtherVertexInEdge1 = verticesByNumber.at(otherVertexInEdge1->Number);
				MeshVertex<2>* coarseOtherVertexInEdge2 = verticesByNumber.at(otherVertexInEdge2->Number);

				DimVector<2> vector1 = Vect<2>(coarseV, coarseOtherVertexInEdge1);
				DimVector<2> vector2 = Vect<2>(coarseV, coarseOtherVertexInEdge2);

				double cosine = vector1.dot(vector2) / (vector1.norm()*vector2.norm());
				if (cosine == -1) // the edges are colinear
				{
					Edge* mergedEdge = new Edge(faceNumber++, coarseOtherVertexInEdge1, coarseOtherVertexInEdge2);
					
					if (face1->IsDomainBoundary && face2->IsDomainBoundary)
						mergedEdge->IsDomainBoundary = true;
					else if (!face1->IsDomainBoundary && !face2->IsDomainBoundary)
						mergedEdge->IsDomainBoundary = false;
					else
						assert(false && "Not sure what to do yet");

					face1->CoarseFace = mergedEdge;
					face2->CoarseFace = mergedEdge;
					mergedEdge->FinerFaces.push_back(face1);
					mergedEdge->FinerFaces.push_back(face2);

					coarseSkeleton->AddFace(mergedEdge);

					coarseOtherVertexInEdge1->Faces.push_back(mergedEdge);
					coarseOtherVertexInEdge2->Faces.push_back(mergedEdge);

					break;
				}
			}
		}

		if (!face1->CoarseFace)
			Utils::Warning("Face " + to_string(face1->Number) + " has not been coarsened.");
	}
	this->FinalizeCoarsening();
}