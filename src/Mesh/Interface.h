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
	Interface(vector<Face<Dim>*> faces) : _faces(faces)
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

	inline vector<Face<Dim>*> Faces()
	{
		return _faces;
	}
	inline map<Vertex*, set<Face<Dim>*>> MapVertexFaces()
	{
		return _mapVertexFaces;
	}
	inline vector<Vertex*> InteriorVertices()
	{
		return _interiorVertices;
	}
	vector<Vertex*> BoundaryVertices()
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