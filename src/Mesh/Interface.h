#pragma once
#include "Face.h"
#include "2D/Edge.h"
using namespace std;

template <int Dim>
class Interface
{
private:
	vector<Face<Dim>*> _faces;
	map<Vertex*, set<Face<Dim>*>> _mapVertexFaces;
	vector<Vertex*> _interiorVertices;
	vector<Vertex*> _boundaryVertices;
	Face<Dim>* _collapsedFace = nullptr;
public:
	Interface(const vector<Face<Dim>*>& faces) : _faces(faces)
	{
		Init();
	}

	void Init()
	{
		if (Dim == 2)
		{
			for (int i = 0; i < _faces.size(); i++)
			{
				Face<Dim>* fi = _faces[i];
				for (Vertex* v : fi->Vertices())
				{
					_mapVertexFaces[v].insert(fi);
					for (int j = i + 1; j < _faces.size(); j++)
					{
						Face<Dim>* fj = _faces[j];
						if (fj->HasVertex(v))
						{
							_interiorVertices.push_back(v);
							_mapVertexFaces[v].insert(fj);
							break;
						}
					}
				}
			}
			
			for (auto const& it : _mapVertexFaces)
			{
				if (it.second.size() == 1)
					_boundaryVertices.push_back(it.first);
				else
					assert(it.second.size() == 2);
			}
		}
	}

	inline vector<Face<Dim>*> Faces() const
	{
		return _faces;
	}

	bool IsPhysicalBoundary()
	{
		if (_faces[0]->IsDomainBoundary)
			return true;
		return _faces[0]->Element1->PhysicalGroupId != _faces[0]->Element2->PhysicalGroupId;
	}

	inline map<Vertex*, set<Face<Dim>*>> MapVertexFaces()
	{
		return _mapVertexFaces;
	}
	inline vector<Vertex*>& InteriorVertices()
	{
		return _interiorVertices;
	}
	vector<Vertex*>& BoundaryVertices()
	{
		return _boundaryVertices;
	}
	inline Face<Dim>* CollapsedFace()
	{
		if (!_collapsedFace)
			CreateCollapsedFace();
		return _collapsedFace;
	}

	bool IsClosingItself()
	{
		return _interiorVertices.size() == _faces.size();
	}

	bool HasHoles()
	{
		return _interiorVertices.size() < _faces.size() - 1;
	}

	FaceCollapsingStatus AnalyzeCollapsing()
	{
		if (_faces.size() < 2)
			return FaceCollapsingStatus::NotEnoughFaces;

		if (HasHoles())
			return FaceCollapsingStatus::InterfaceHasHoles;

		Element<Dim>* e1 = _faces[0]->Element1;
		Element<Dim>* e2 = _faces[0]->Element2;

		if (e1->Faces.size() - _faces.size() + 1 <= Dim)
			return FaceCollapsingStatus::ElementFullDegeneration;
		if (e2 && e2->Faces.size() - _faces.size() + 1 <= Dim)
			return FaceCollapsingStatus::ElementFullDegeneration;

		assert(_boundaryVertices.size() == 2);

		Element<Dim>* smallElem = !e2 || e1->Faces.size() < e2->Faces.size() ? e1 : e2;
		Element<Dim>* bigElem = smallElem == e1 ? e2 : e1;

		// Analyze the small element first because it has more chance to degenerate
		FaceCollapsingStatus status = AnalyzeCollapsing(smallElem);
		if (status != FaceCollapsingStatus::Ok)
			return status;
		if (e2)
		{
			status = AnalyzeCollapsing(bigElem);
			if (status != FaceCollapsingStatus::Ok)
				return status;

			if (!e1->IsConvex() && e1->ConvexHullEmbeds(e2))
			{
				//PlotMatlab();
				return FaceCollapsingStatus::OneElementEmbeddedInConvexHullOfTheOther;
			}
			if (!e2->IsConvex() && e2->ConvexHullEmbeds(e1))
			{
				//PlotMatlab();
				return FaceCollapsingStatus::OneElementEmbeddedInConvexHullOfTheOther;
			}
		}

		return FaceCollapsingStatus::Ok;
	}

private:
	FaceCollapsingStatus AnalyzeCollapsing(Element<Dim>* element)
	{
		if (!_collapsedFace)
			CreateCollapsedFace();

		for (Face<Dim>* f : element->Faces)
		{
			// If a face which is not part of the merged faces...
			if (!f->IsIn(_faces))
			{
				// ... is included in the collaped face
				if (_collapsedFace->Contains(f->Center()))
				{
					// then the element partially degenerates.
					delete _collapsedFace;
					return FaceCollapsingStatus::ElementPartialDegeneration;
				}

				// ... and is not connected to the collapsed face by its extremities
				if (!f->HasCommonVerticesWith(_collapsedFace))
				{
					/*for (Vertex* v : f->Vertices())
					{
						// ... but has a vertex contained in the collapsed face
						if (collapsedFace->Contains(*v))
							// then the polygon is split in two parts connected by a single vertex.
							return true;
					}*/
					// ... and crosses the collapsed face
					if (f->IntersectsWith(_collapsedFace))
					{
						delete _collapsedFace;
						return FaceCollapsingStatus::CrossedPolygon;
					}
				}
			}
		}
		return FaceCollapsingStatus::Ok;
	}

public:
	list<set<Face<Dim>*>> CoplanarSubsets() const
	{
		assert(Dim == 2);

		list<set<Face<Dim>*>> subsets;

		for (Vertex* v : _interiorVertices)
		{
			set<Face<Dim>*> vFaces = _mapVertexFaces.at(v);
			auto it = vFaces.cbegin();
			Face<Dim>* f1 = *it; it++;
			Face<Dim>* f2 = *it;
			vector<Vertex*> f1Vertices = f1->Vertices();
			vector<Vertex*> f2Vertices = f2->Vertices();

			DimVector<2> v1 = Vect<2>(f1Vertices[0], f1Vertices[1]);
			DimVector<2> v2 = Vect<2>(f2Vertices[0], f2Vertices[1]);
			if (AreCollinear(v1, v2))
				subsets.push_back({ f1, f2 });
		}

		// Iterate over the pairs of collinear faces
		for (auto itSubset1 = subsets.begin(); itSubset1 != subsets.end(); itSubset1++)
		{
			set<Face<Dim>*>& subset1 = *itSubset1;

			// Iterate over the current subsets to see if it can be put into one
			for (auto itSubset2 = subsets.begin(); itSubset2 != subsets.end(); )
			{
				if (itSubset1 == itSubset2)
				{
					itSubset2++;
					continue;
				}

				set<Face<Dim>*> subset2 = *itSubset2;
				bool areDisjoint = true;
				for (auto it = subset2.begin(); it != subset2.end(); )
				{
					if (subset1.find(*it) != subset1.end())
					{
						areDisjoint = false;
						break;
					}
					it++;
				}

				if (!areDisjoint)
				{
					subset1.insert(subset2.begin(), subset2.end()); // merge
					subsets.erase(itSubset2);
					itSubset2 = subsets.begin();
				}
				else
					itSubset2++;
			}
		}

		return subsets;
	}

private:
	void PlotMatlab()
	{
		Element<Dim>* e1 = _faces[0]->Element1;
		Element<Dim>* e2 = _faces[0]->Element2;

		cout << "%--------------------- Elem1:" << endl;
		e1->ExportToMatlab();
		if (e2)
		{
			cout << "%--------------------- Elem2:" << endl;
			e2->ExportToMatlab("b");
		}
		cout << "%--------------------- MergedFace:" << endl;
		_collapsedFace->ExportToMatlab("m");
	}

	void CreateCollapsedFace() { assert(false); }
};

template<>
void Interface<2>::CreateCollapsedFace()
{
	_collapsedFace = new Edge(0, _boundaryVertices[0], _boundaryVertices[1]);
}