#pragma once
#include <set>
#include <list>
#include <algorithm>
#include "Mesh.h"
#include "2D/Polygon.h"
#include "AgglomerateElement.h"
using namespace std;

template <int Dim>
class PolyhedralMesh : public Mesh<Dim>
{
private:
	double _regularity = 0;
public:
	const int AgglomerateSize = 4*Dim;

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
		else if (strategy == CoarseningStrategy::AgglomerationCoarseningByClosestCenter || strategy == CoarseningStrategy::AgglomerationCoarseningByLargestInterface)
			CoarsenByAgglomerationByPairs(strategy);
		else if (strategy == CoarseningStrategy::AgglomerationCoarseningBySeedPoints)
			CoarsenByAgglomerationBySeedPoints();
		else if (strategy == CoarseningStrategy::FaceCoarsening)
			FaceCoarsening();
		else
			Mesh<Dim>::CoarsenMesh(strategy);

		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
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

		BigNumber elementToAnalyze = 999999999999;

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

			if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze)
			{
				coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
				AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
				AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
				shape->ExportSubShapesToMatlab();
			}

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

				if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze)
				{
					coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
					AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
					AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
					shape->ExportSubShapesToMatlab();
				}

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
			if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze)
			{
				coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
				AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
				AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
				shape->ExportSubShapesToMatlab();
			}

			coarseMesh->AgglomerateFineFaces(f1, f2);

			if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze)
			{
				coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
				AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
				AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
				shape->ExportSubShapesToMatlab();
			}

			coarseMesh->CollapseInterfacesMadeOfMultipleFaces(coarseElement1);
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

			if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze && (e == coarseMesh->Elements[elementToAnalyze] || smallestNeighbour == coarseMesh->Elements[elementToAnalyze]))
			{
				coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
				AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
				AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
				shape->ExportSubShapesToMatlab();
			}
		}


		//------------------------------------------------------------------------------//
		// We merge the faces which are part of the same interface between two elements //
		//------------------------------------------------------------------------------//
		for (Element<Dim>* ce : coarseMesh->Elements)
			coarseMesh->CollapseInterfacesMadeOfMultipleFaces(ce);

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

	void CollapseInterfacesMadeOfMultipleFaces(Element<Dim>* e)
	{
		// avoid the destruction of the element
		if (e->Faces.size() == 3)
			return;

		for (Face<Dim>* f : e->Faces)
		{
			if (f->IsDomainBoundary)
				continue;

			Element<Dim>* n = f->GetNeighbour(e);
			if (!n) // no neighbour is affected yet
				continue;

			// avoid the destruction of the neighbour
			if (n->Faces.size() == 3)
				continue;

			vector<Face<Dim>*> faces = e->InterfaceWith(n);
			if (faces.size() > 1)
			{
				/*if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze && (ce == coarseMesh->Elements[elementToAnalyze] || n == coarseMesh->Elements[elementToAnalyze]))
					//if (this->FineMesh && this->FineMesh->FineMesh && faces[0]->Number == 289 && faces[1]->Number == 292)
				{
					coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse3.dat");
					AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
					AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
					shape->ExportSubShapesToMatlab();
				}*/
				/*if (this->FineMesh && this->FineMesh->FineMesh && n->Number == 136)
				{
					this->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse3.dat");
					AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(this->Elements[136]);
					AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
					shape->ExportSubShapesToMatlab();
				}*/
				this->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
				Agglomerate(faces);

				/*if (this->FineMesh && coarseMesh->Elements.size() > elementToAnalyze && (ce == coarseMesh->Elements[elementToAnalyze] || n == coarseMesh->Elements[elementToAnalyze]))
				{
					coarseMesh->ExportFacesToMatlab("/mnt/c/Users/pierr/Desktop/Matrices/coarse2.dat");
					AgglomerateElement<Dim>* agglo = dynamic_cast<AgglomerateElement<Dim>*>(coarseMesh->Elements[elementToAnalyze]);
					AgglomerateShape<Dim>* shape = dynamic_cast<AgglomerateShape<Dim>*>(agglo->Shape());
					shape->ExportSubShapesToMatlab();
				}*/

				break;
			}
		}
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
		/*if (this->Elements.size() <= AgglomerateSize)
		{
			cout << "Error: impossible to build coarse mesh. Only " << this->Elements.size() << " element(s) left." << endl;
		}
		else*/
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

			// Find the closest neighbour of e still available for aggregation
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
				double distance = Vect<2>(e->Center(), neighbour->Center()).norm();
				double interfaceMeasure = 0;
				if (strategy == CoarseningStrategy::AgglomerationCoarseningByLargestInterface)
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
					else if (strategy == CoarseningStrategy::AgglomerationCoarseningByLargestInterface && interfaceMeasure > largestInterface)
					{
						neighbourForAggreg = neighbour;
						largestInterface = interfaceMeasure;
					}
				}
			}

			if (neighbourForAggreg)
			{
				// Agglomeration
				Element<Dim>* macroElement = coarseMesh->AgglomerateFineElements({ e, neighbourForAggreg });
				coarseMesh->CollapseInterfacesMadeOfMultipleFaces(macroElement);
			}
			else
			{
				// If all the neighbours have already been aggregated, then we aggregate with the closest macroElement.
				AgglomerateElement<Dim>* coarseNeighbourForAggreg = nullptr;
				for (Face<Dim>* f : e->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					AgglomerateElement<Dim>* macroNeighbour = dynamic_cast<AgglomerateElement<Dim>*>(f->GetNeighbour(e)->CoarserElement);
					double distance = Vect<2>(e->Center(), macroNeighbour->Center()).norm();
					if (!coarseNeighbourForAggreg || distance < closestDistance)
					{
						coarseNeighbourForAggreg = macroNeighbour;
						closestDistance = distance;
					}
				}

				if (!coarseNeighbourForAggreg)
					Utils::FatalError("Element cannot be aggregated. Weird...");

				// Agglomeration
				coarseMesh->AgglomerateFineElementToCoarse(e, coarseNeighbourForAggreg);
				coarseMesh->CollapseInterfacesMadeOfMultipleFaces(coarseNeighbourForAggreg);
			}
		}

		for (Element<Dim>* ce : coarseMesh->Elements)
			coarseMesh->CollapseInterfacesMadeOfMultipleFaces(ce);
	}



	//--------------------------------------//
	//                                      //
	//        General useful methods        //
	//                                      //
	//--------------------------------------//


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

		// Collect the faces interfacing the elements
		vector<Face<Dim>*> facesToRemove;

		// Initialisation of the coarse macro-element
		AgglomerateElement<Dim>* coarseElement = new AgglomerateElement<Dim>(this->Elements.size(), fineElements[0]);
		coarseElement->Faces = fineElements[0]->Faces;
		fineElements[0]->CoarserElement = coarseElement;

		for (int i = 0; i < fineElements.size(); i++)
		{
			Element<Dim>* e1 = fineElements[i];
			for (Face<Dim>* f : e1->Faces)
			{
				Element<Dim>* neighbour = f->GetNeighbour(e1);
				bool neighbourMustBeAgglomerated = false;
				for (int j = i + 1; j < fineElements.size(); j++)
				{
					if (fineElements[j] == neighbour)
					{
						neighbourMustBeAgglomerated = true;
						break;
					}
				}

				if (neighbourMustBeAgglomerated)
				{
					facesToRemove.push_back(f);
					if (!neighbour->CoarserElement)
					{
						// Agglomeration
						coarseElement->Add(neighbour);
						neighbour->CoarserElement = coarseElement;
						coarseElement->Faces = Utils::Join(neighbour->NonCommonFacesWith(coarseElement), coarseElement->NonCommonFacesWith(neighbour));
					}
				}
			}
		}

		for (Element<Dim>* e : fineElements)
		{
			e->CoarserElement = coarseElement;
			coarseElement->FinerElements.push_back(e);
		}

		for (Face<Dim>* f : facesToRemove)
		{
			f->IsRemovedOnCoarserGrid = true;
			coarseElement->FinerFacesRemoved.push_back(f);
		}

		// Add coarse element to the list
		coarseElement->Number = this->Elements.size();
		this->Elements.push_back(coarseElement);

		// Kept faces are cloned for the coarse mesh and linked to their clones and the coarse element.
		vector<Face<Dim>*> facesToClone = coarseElement->Faces;
		coarseElement->Faces.clear();
		for (Face<Dim>* f : facesToClone)
			this->CloneAndAddFace(f, coarseElement);

		coarseElement->Init();
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
				this->RemoveFace(f->CoarseFace);
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
			this->RemoveFace(ff->CoarseFace);
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

	void AgglomerateFineFaces(Face<Dim>* f1, Face<Dim>* f2)
	{
		// f1 and f2 are one level higher than 'this'
		Agglomerate({ f1->CoarseFace, f2->CoarseFace });
	}

	void Agglomerate(vector<Face<Dim>*> faces)
	{
		assert(faces.size() > 1);
		// 'faces' must be at the same level as 'this'

		// Identification of the vertices to remove (those which belong to 2 faces)
		vector<Vertex*> verticesToRemove;
		map<Vertex*, set<Face<Dim>*>> mapVertexFaces;
		for (int i = 0; i < faces.size(); i++)
		{
			Face<Dim>* fi = faces[i];
			for (Vertex* v : fi->Vertices())
			{
				mapVertexFaces[v].insert(fi);
				for (int j = i + 1; j < faces.size(); j++)
				{
					Face<Dim>* fj = faces[j];
					if (fj->HasVertex(v))
					{
						verticesToRemove.push_back(v);
						mapVertexFaces[v].insert(fj);
						break;
					}
				}
			}
		}
		if (verticesToRemove.size() != faces.size() - 1)
		{
			Utils::Warning("Attempt to perform an operation which would result in an element rounding another one. Refused.");
			return;
		}
		//assert(verticesToRemove.size() == faces.size() - 1);

		// The two vertices to keep
		Vertex* v1 = nullptr;
		Vertex* v2 = nullptr;

		for (auto const& it : mapVertexFaces)
		{
			if (it.second.size() == 1)
			{
				if (!v1)
					v1 = it.first;
				else if (!v2)
				{
					v2 = it.first;
					break;
				}
			}
			else
				assert(it.second.size() == 2);
		}

		Edge* mergedEdge = new Edge(0, v1, v2);
		ReplaceFaces(faces, dynamic_cast<Face<Dim>*>(mergedEdge));
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
			this->RemoveFace(f);
			delete f;
		}

		// Add mergedFace to the list
		this->AddFace(mergedFace);
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
			this->AddFace(copy);
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
	Polygon* macroElement = new Polygon(elementNumber, vertices);
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
	Polygon* macroElement = new Polygon(0, e1, e2, facesToRemove);
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