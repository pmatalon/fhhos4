#pragma once
#include <set>
#include <list>
#include <algorithm>
#include "Mesh.h"
#include "2D/PolygonalElement.h"
#include "AgglomerateElement.h"
#include "Agglo.h"
using namespace std;

template <int Dim>
class PolyhedralMesh : public Mesh<Dim>
{
private:
	double _regularity = 0;
public:
	PolyhedralMesh() : Mesh<Dim>()
	{}

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

	virtual double H() override
	{
		assert(false);
	}

	virtual double Regularity() override
	{
		if (_regularity == 0)
		{
			for (Element<Dim>* e : this->Elements)
			{
				if (e->Regularity() < _regularity)
					_regularity = e->Regularity();
			}
		}
		return _regularity;
	}

	virtual void CoarsenMesh(CoarseningStrategy strategy) override
	{
		if (this->CoarseMesh)
			return;
		if (strategy == CoarseningStrategy::SplittingRefinement || strategy == CoarseningStrategy::BeyRefinement)
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
			Mesh<Dim>::CoarsenMesh(strategy);

		this->CoarseMesh->FillBoundaryAndInteriorFaceLists();

		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
	}

	virtual void Init()
	{
		ElementParallelLoop<Dim> parallelLoop(this->Elements);
		parallelLoop.Execute([](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				// Make Init() available for all element classes
				PolygonalShape* p = dynamic_cast<PolygonalShape*>(e->Shape());
				if (p)
				{
					p->ComputeTriangulation();
					p->ComputeBoundingBox();
				}
			});
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
		this->CoarseMesh = coarseMesh;
		coarseMesh->FineMesh = this;

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

			/*if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze)
			{
				coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
				AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
				AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
				shape->ExportSubShapesToMatlab();
			}*/

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

				/*if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze)
				{
					coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
					AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
					AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
					shape->ExportSubShapesToMatlab();
				}*/

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
			/*if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze)
			{
				coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
				AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
				AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
				shape->ExportSubShapesToMatlab();
			}*/

			FaceCollapsingStatus status = coarseMesh->AgglomerateFineFaces(f1, f2);
			assert(status == FaceCollapsingStatus::Ok);

			/*if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze)
			{
				coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
				AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
				AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
				shape->ExportSubShapesToMatlab();
			}*/

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
				coarseMesh->AgglomerateFineElementToCoarse(e, dynamic_cast<AgglomerateElement<Dim>*>(smallestNeighbour));

			/*if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze && (e == coarseMesh->Elements[elementToAnalyze] || smallestNeighbour == coarseMesh->Elements[elementToAnalyze]))
			{
				coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
				AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
				AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
				shape->ExportSubShapesToMatlab();
			}*/
		}


		//------------------------------------------------------------------------------//
		// We merge the faces which are part of the same interface between two elements //
		//------------------------------------------------------------------------------//
		for (Element<Dim>* ce : coarseMesh->Elements)
			coarseMesh->TryCollapseInterfacesMadeOfMultipleFaces(ce);

		/*if (this->FineMesh && coarseMesh->Elements.size() > 8)
		{
			coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
			AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[8]);
			AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
			shape->ExportSubShapesToMatlab();
		}*/

		// Set the coarse vertex list
		set<Vertex*> vertices;
		for (auto f : coarseMesh->Faces)
		{
			for (auto v : f->Vertices())
				vertices.insert(v);
		}
		coarseMesh->Vertices = vector<Vertex*>(vertices.begin(), vertices.end());
	}

	//-------------------------------------------------------------------------------------------------------------------------//
	// Browse the neighbours: if the interface with one of them is made of multiple faces, they're collapsed into a single one //
	//-------------------------------------------------------------------------------------------------------------------------//
	void TryCollapseInterfacesMadeOfMultipleFaces(Element<Dim>* e)
	{
		for (Element<Dim>* n : e->Neighbours())
		{
			FaceCollapsingStatus status = TryCollapseInterfaceBetween(e, n);
			if (status == FaceCollapsingStatus::ElementFullDegeneration || status == FaceCollapsingStatus::OneElementEmbeddedInConvexHullOfTheOther)
			{
				// We can't collpase the faces, otherwise one element will degenerate. Agglomerating the elements instead.
				Element<Dim>* newE = Agglomerate(e, n);
				TryCollapseInterfacesMadeOfMultipleFaces(newE);
				return;
			}
			// Any other failing status, abort the collapsing.
		}

		/*if (e->IsOnBoundary())
		{
			// TODO agglomerate (at least) collinear faces
		}*/
	}

	FaceCollapsingStatus TryCollapseInterfaceBetween(Element<Dim>* e1, Element<Dim>* e2)
	{
		vector<Face<Dim>*> interfaceFaces = e1->InterfaceWith(e2);
		if (interfaceFaces.size() > 1)
			return TryCollapse(interfaceFaces);
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
		this->CoarseMesh = coarseMesh;
		coarseMesh->FineMesh = this;

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
		this->CoarseMesh = coarseMesh;
		coarseMesh->FineMesh = this;

		for (Element<Dim>* currentElem : this->Elements)
		{
			if (currentElem->CoarserElement)
				continue;

			vector<Element<Dim>*> availableNeighbours = AvailableFaceNeighbours(currentElem);

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

				Element<Dim>* coarseElement = coarseMesh->AgglomerateFineElementToCoarse(currentElem, coarseNeighbourForAggreg);
				coarseMesh->TryCollapseInterfacesMadeOfMultipleFaces(coarseElement);
			}
		}

		list<Face<Dim>*> uncoarsenedFaces;
		for (Face<Dim>* f : this->Faces)
		{
			if (!f->IsDomainBoundary && !f->IsRemovedOnCoarserGrid && !f->HasBeenCoarsened())
				uncoarsenedFaces.push_back(f);
		}

		/*cout << uncoarsenedFaces.size() << " faces not removed or coarsened remain." << endl;

		typename list<Face<Dim>*>::iterator it = uncoarsenedFaces.begin();
		while (it != uncoarsenedFaces.end())
		{
			Face<Dim>* f = *it;
			if (f->IsRemovedOnCoarserGrid || f->HasBeenCoarsened())
				uncoarsenedFaces.erase(it);
			else
			{
				TryCollapseInterfaceBetween(f->Element1, f->Element2);
				if (interfaceHasBeenCollapsed)
					uncoarsenedFaces.erase(it);
			}
			it++;
		}*/

		cout << uncoarsenedFaces.size() << " faces not removed or coarsened remain (out of " << this->Faces.size() << " faces)." << endl;

		coarseMesh->Init();
	}

	//-----------------------------------------------//
	//                                               //
	//        Coarsening by vertex neighbours        //
	//                                               //
	//-----------------------------------------------//

	void CoarsenByAgglomerationByVertexNeighbours()
	{
		PolyhedralMesh<Dim>* coarseMesh = new PolyhedralMesh<Dim>();
		this->CoarseMesh = coarseMesh;
		coarseMesh->FineMesh = this;

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
		coarseMesh->Init();
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
				for (Vertex* v2 : e->Shape()->Vertices())
					vertexElements[v2].clear();
			}
		}

		this->CoarseMesh = coarseMesh;
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
			for (Vertex* v : e->Shape()->Vertices())
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
	Element<Dim>* AgglomerateByVertexRemoval(vector<Element<Dim>*> fineElements, Vertex* removedVertex)
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
						for (Vertex* v : f->Shape()->Vertices())
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
		macroElement->Number = this->Elements.size();
		this->Elements.push_back(macroElement);

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
		coarseMesh->Init();
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
		this->CoarseMesh = coarseMesh;
		coarseMesh->FineMesh = this;

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

	vector<Element<Dim>*> AvailableFaceNeighbours(Element<Dim>* elem)
	{
		vector<Element<Dim>*> availableNeighbours;
		for (Element<Dim>* n : elem->Neighbours())
		{
			if (!n->CoarserElement)
				availableNeighbours.push_back(n);
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
		Element<Dim>* coarseNeighbourForAggreg = nullptr;
		double closestDistance = -1;
		double largestInterface = -1;
		for (Face<Dim>* f : e->Faces)
		{
			if (f->IsDomainBoundary)
				continue;

			Element<Dim>* macroNeighbour = f->GetNeighbour(e)->CoarserElement;

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


	//--------------------------------------------------------------------//
	// Creates the macro-element obtained by agglomeration of 2 elements, //
	// adds it to the mesh, and return it.                                //
	//--------------------------------------------------------------------//
	Element<Dim>* AgglomerateFineElements(Element<Dim>* e1, Element<Dim>* e2)
	{
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
		coarseElement->Number = this->Elements.size();
		this->Elements.push_back(coarseElement);

		// Kept faces are cloned for the coarse mesh and linked to their clones and the coarse element.
		for (Face<Dim>* f : facesToClone)
			this->CloneAndAddFace(f, coarseElement);

		assert(coarseElement->Faces.size() > 2);
		if (Dim == 2)
			assert(coarseElement->Faces.size() == coarseElement->Vertices().size());

		return coarseElement;
	}

	Element<Dim>* AgglomerateFineElements(vector<Element<Dim>*> fineElements)
	{
		//assert(fineElements.size() > 1 || fineElements[0]->IsOnBoundary());
		for (Element<Dim>* e : fineElements)
			assert(!e->CoarserElement);

		Agglo<Dim> agglo(fineElements);
		Element<Dim>* coarseElement = CreatePolyhedron(agglo.Vertices());

		for (Element<Dim>* e : fineElements)
		{
			e->CoarserElement = coarseElement;
			coarseElement->FinerElements.push_back(e);
		}

		for (Face<Dim>* f : agglo.RemovedFaces())
		{
			f->IsRemovedOnCoarserGrid = true;
			coarseElement->FinerFacesRemoved.push_back(f);
		}

		// Add coarse element to the list
		coarseElement->Number = this->Elements.size();
		this->Elements.push_back(coarseElement);

		// Kept faces are cloned for the coarse mesh and linked to their clones and the coarse element.
		for (Face<Dim>* f : agglo.Faces)
			this->CloneAndAddFace(f, coarseElement);

		assert(coarseElement->Vertices().size() > 2);

		assert(coarseElement->Faces.size() > 2);
		if (Dim == 2)
			assert(coarseElement->Faces.size() == coarseElement->Vertices().size());

		return coarseElement;
	}

	void AgglomerateFineElementToCoarse(Element<Dim>* fineElement, AgglomerateElement<Dim>* coarseElement)
	{
		assert(!fineElement->CoarserElement);

		// Add the fine element to the agglomerate
		coarseElement->Add(fineElement);
		fineElement->CoarserElement = coarseElement;
		coarseElement->FinerElements.push_back(fineElement);

		// Collect the faces interfacing the elements
		for (Face<Dim>* f : fineElement->Faces)
		{
			// Face to remove
			if (!f->IsDomainBoundary && f->GetNeighbour(fineElement)->CoarserElement == coarseElement)
			{
				f->IsRemovedOnCoarserGrid = true;
				coarseElement->FinerFacesRemoved.push_back(f);
				coarseElement->RemoveFace(f->CoarseFace);
				this->RemoveFace(f->CoarseFace, false);
				delete f->CoarseFace;
				f->CoarseFace = nullptr;
			}
			else // Kept faces are cloned for the coarse mesh and linked to their clones and the coarse element.
				this->CloneAndAddFace(f, coarseElement);
		}

		// Recompute the shape of the agglomerate
		coarseElement->Init();
		assert(coarseElement->Vertices().size() > 2);
		assert(coarseElement->Faces.size() > 2);
		if (Dim == 2)
			assert(coarseElement->Faces.size() == coarseElement->Vertices().size());
	}


	Element<Dim>* AgglomerateFineElementToCoarse(Element<Dim>* fineElement, Element<Dim>* coarseElement)
	{
		// Collect the fine faces interfacing these two elements //
		vector<Face<Dim>*> facesToRemove;
		vector<Face<Dim>*> facesToClone;
		vector<Face<Dim>*> coarseFacesToKeep;
		for (Face<Dim>* ff : fineElement->Faces)
		{
			if (!ff->IsDomainBoundary && ff->GetNeighbour(fineElement)->CoarserElement == coarseElement)
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
		
		// Creation of the polygonal macro-element
		Element<Dim>* newCoarseElement = CreateMacroElement(fineElement, coarseElement, facesToRemove);

		// Links between the macro element and the fine elements
		for (Element<Dim>* fe : coarseElement->FinerElements)
		{
			newCoarseElement->FinerElements.push_back(fe);
			fe->CoarserElement = newCoarseElement;
		}
		newCoarseElement->FinerElements.push_back(fineElement);
		fineElement->CoarserElement = newCoarseElement;

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

			// Remove f->CoarseFace from this->Faces
			this->RemoveFace(ff->CoarseFace, false);
			delete ff->CoarseFace;
			ff->CoarseFace = nullptr;
		}

		// Delete the old coarse element and replace it with the new one
		newCoarseElement->Number = coarseElement->Number;
		delete coarseElement;
		this->Elements[newCoarseElement->Number] = newCoarseElement;

		// Kept faces are cloned for the coarse mesh and linked to their clones and the coarse element.
		for (Face<Dim>* ff : facesToClone)
			this->CloneAndAddFace(ff, newCoarseElement);

		assert(newCoarseElement->Faces.size() > 2);
		if (Dim == 2)
			assert(newCoarseElement->Faces.size() == newCoarseElement->Vertices().size());

		return newCoarseElement;
	}

	Element<Dim>* Agglomerate(Element<Dim>* e1, Element<Dim>* e2)
	{
		vector<Element<Dim>*> v{ e1, e2 };
		return Agglomerate(v);
	}

	Element<Dim>* Agglomerate(vector<Element<Dim>*> elements)
	{
		Agglo<Dim> agglo(elements);
		Element<Dim>* newElement = CreatePolyhedron(agglo.Vertices());

		for (Element<Dim>* e : elements)
		{
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
			this->RemoveFace(f, false);
		}

		newElement->Faces = agglo.Faces;

		// Replace the first element with the new one
		newElement->Number = elements[0]->Number;
		this->Elements[newElement->Number] = newElement;
		delete elements[0];
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

	FaceCollapsingStatus TryCollapse(vector<Face<Dim>*> faces)
	{
		assert(faces.size() > 1);
		// 'faces' must be at the same level as 'this'

		Interface<Dim> interf(faces);

		FaceCollapsingStatus status = interf.AnalyzeCollapsing();
		if (status != FaceCollapsingStatus::Ok)
			return status;

		ReplaceFaces(interf.Faces(), interf.CollapsedFace());
		return FaceCollapsingStatus::Ok;
	}

	void ReplaceFaces(vector<Face<Dim>*> faces, Face<Dim>* mergedFace)
	{
		// The faces and mergedFace must be at the same level.

		mergedFace->IsDomainBoundary = faces[0]->IsDomainBoundary;

		// Replace in the elements
		Element<Dim>* elem1 = faces[0]->Element1;
		Element<Dim>* elem2 = faces[0]->Element2;
		mergedFace->Element1 = elem1;
		mergedFace->Element2 = elem2;

		elem1->ReplaceFaces(faces, mergedFace);
		if (elem2)
			elem2->ReplaceFaces(faces, mergedFace);

		// Links
		for (Face<Dim>* f : faces)
		{
			for (Face<Dim>* ff : f->FinerFaces)
			{
				mergedFace->FinerFaces.push_back(ff);
				ff->CoarseFace = mergedFace;
			}
		}

		// Remove the faces from the mesh
		for (Face<Dim>* f : faces)
		{
			this->RemoveFace(f, false);
			delete f;
		}

		// Add mergedFace to the list
		this->AddFace(mergedFace, false);
	}

	void CloneAndAddFace(Face<Dim>* f, Element<Dim>* macroElement)
	{
		if (!f->CoarseFace) // this face has not already been created by a previous aggregation
		{
			// Clone face
			BigNumber faceNumber = this->Faces.size();
			Face<Dim>* copy = f->CreateSameGeometricFace(faceNumber, macroElement);
			copy->IsDomainBoundary = f->IsDomainBoundary;

			// Associate
			f->CoarseFace = copy;
			copy->FinerFaces.push_back(f);

			// Add copied face to the list of coarse faces
			this->AddFace(copy, false);
		}
		else
		{
			assert(f->CoarseFace->Element2 == nullptr);
			f->CoarseFace->Element2 = macroElement;

			assert(f->CoarseFace->Element1 != f->CoarseFace->Element2);
		}
		macroElement->AddFace(f->CoarseFace);
	}

	//-----------------------------------------------//
	// Dim-specific functions (implementation below) //
	//-----------------------------------------------//
	Element<Dim>* CreatePolyhedron(vector<Vertex*> vertices) { return nullptr; }
	Element<Dim>* CreateMacroElement(Element<Dim>* e1, Element<Dim>* e2, vector<Face<Dim>*> facesToRemove) { return nullptr; }
	Face<Dim>* CreateMacroFace(Face<Dim>* f1, Face<Dim>* f2, Vertex* vertexToRemove) { return nullptr; }
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
Element<2>* PolyhedralMesh<2>::CreateMacroElement(Element<2>* e1, Element<2>* e2, vector<Face<2>*> facesToRemove)
{
	PolygonalElement* macroElement = new PolygonalElement(0, e1, e2, facesToRemove, false);
	return macroElement;
}

template<>
Element<3>* PolyhedralMesh<3>::CreateMacroElement(Element<3>* e1, Element<3>* e2, vector<Face<3>*> facesToRemove)
{
	Utils::FatalError("Not implemented in 3D.");
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
}

template <>
void PolyhedralMesh<2>::FaceCoarsening()
{
	PolyhedralMesh<2>* coarseSkeleton = new PolyhedralMesh<2>();
	coarseSkeleton->ComesFrom.CS = CoarseningStrategy::FaceCoarsening;
	this->CoarseMesh = coarseSkeleton;
	coarseSkeleton->FineMesh = this;

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
}