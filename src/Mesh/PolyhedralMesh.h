#pragma once
#include <set>
#include <list>
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

	virtual void CoarsenMesh(CoarseningStrategy strategy) override
	{
		if (strategy == CoarseningStrategy::Standard)
			CoarsenByAgglomerationAndMergeColinearFaces();
		else if (strategy == CoarseningStrategy::Agglomeration)
			CoarsenByAgglomerationAndKeepFineFaces();
		else
			assert(false && "Coarsening strategy not implemented!");
		this->CoarseMesh->SetDiffusionCoefficient(this->_diffusionPartition);
		this->CoarseMesh->SetBoundaryConditions(this->_boundaryConditions);
	}

private:

	void CoarsenByAgglomerationAndMergeColinearFaces()
	{
	}

	void CoarsenByAgglomerationAndKeepFineFaces()
	{
		if (this->Elements.size() <= AgglomerateSize)
		{
			cout << "Error: impossible to build coarse mesh. Only " << this->Elements.size() << " element(s) left." << endl;
		}
		else
		{
			PolyhedralMesh<Dim>* coarseMesh = new PolyhedralMesh<Dim>();

			// TODO: iterate on the interior vertices instead of the elements
			list<Element<Dim>*> remainingElements;// (this->Elements);
			while (!remainingElements.empty())
			{
				Element<Dim>* e = remainingElements.front();
				remainingElements.pop_front();
				if (!e->CoarserElement) // not already aggregated
				{
					vector<Element<Dim>*> elementsToAggregate;
					elementsToAggregate.push_back(e);
					for (Face<Dim>* f : e->Faces)
					{
						if (!f->IsDomainBoundary)
						{

						}
					}
					coarseMesh->Agglomerate(elementsToAggregate);
				}
			}

			this->CoarseMesh = coarseMesh;
		}
	}

private:

	Element<Dim>* Agglomerate(vector<Element<Dim>*> fineElements)
	{
		// The interior faces of the future macro-element are flagged.
		for (Element<2>* e1 : fineElements)
		{
			for (Element<2>* e2 : fineElements)
			{
				if (e1 != e2)
				{
					Face<2>* interface = e1->InterfaceWith(e2);
					interface->IsRemovedOnCoarserGrid = true;
				}
			}
		}

		// Separation of the faces we keep and the faces we remove
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
				Face<Dim>* copy = f->CreateSameGeometricFace(0, macroElement);
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

	Element<Dim>* CreatePolyhedron(vector<Vertex*> vertices)
	{
		if (Dim == 2)
		{
			Polygon* macroElement = new Polygon(0, vertices);
			return macroElement;
		}
		assert(false && "Not yet implemented");
	}
};