#pragma once
#include <set>
#include <list>
#include <algorithm>
#include "Mesh.h"
#include "2D/Polygon.h"
using namespace std;

template <int Dim>
class PolyhedralMesh : public Mesh<Dim>
{
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
		assert(false);
	}

	virtual void CoarsenMesh(CoarseningStrategy strategy) override
	{
		if (this->CoarseMesh)
			return;
		if (strategy == CoarseningStrategy::SplittingRefinement || strategy == CoarseningStrategy::BeyRefinement)
			return;
		if (strategy == CoarseningStrategy::StandardCoarsening)
			CoarsenByAgglomerationAndMergeColinearFaces();
		else if (strategy == CoarseningStrategy::AgglomerationCoarsening)
			CoarsenByAgglomerationAndKeepFineFaces();
		else if (strategy == CoarseningStrategy::FaceCoarsening)
			FaceCoarsening();
		else
			Mesh<Dim>::CoarsenMesh(strategy);

		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
	}

private:

	void CoarsenByAgglomerationAndMergeColinearFaces()
	{
		assert(false);
	}

	void CoarsenByAgglomerationAndKeepFineFaces()
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

				coarseMesh->Agglomerate(elementsAroundV, v);

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

private:
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

	Element<Dim>* Agglomerate(vector<Element<Dim>*> fineElements, Vertex* removedVertex)
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
					Face<Dim>* interface = e1->InterfaceWith(e2);
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

		// Kept faces are cloned for the coarse mesh and linked to their clones.
		for (Face<Dim>* f : keptFaces)
		{
			if (!f->CoarseFace) // this face has not been already created by a previous aggregation
			{
				BigNumber faceNumber = this->Faces.size();
				Face<Dim>* copy = f->CreateSameGeometricFace(faceNumber, macroElement);
				copy->IsDomainBoundary = f->IsDomainBoundary;
				f->CoarseFace = copy;
				copy->FinerFaces.push_back(f);
				// Add copied face to the list of coarse faces
				this->AddFace(copy);
			}
			else
			{
				assert(f->CoarseFace->Element2 == nullptr);
				f->CoarseFace->Element2 = macroElement;
			}
			macroElement->AddFace(f->CoarseFace);
		}

		for (Face<Dim>* f : removedFaces)
			macroElement->FinerFacesRemoved.push_back(f);

		// Add macro-element to the list
		this->Elements.push_back(macroElement);

		return macroElement;
	}

	// Dim-specific functions
	Element<Dim>* CreatePolyhedron(vector<Vertex*> vertices) { return nullptr; }
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