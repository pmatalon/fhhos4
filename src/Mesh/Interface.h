#pragma once
#include "Face.h"
using namespace std;

template <int Dim>
class Interface
{
private:
	vector<Face<Dim>*> _faces;
	map<Vertex*, set<Face<Dim>*>> _mapVertexFaces;
	vector<Vertex*> _interiorVertices;
	vector<Vertex*> _boundaryVertices;
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

	bool HasHoles()
	{
		return _interiorVertices.size() != _faces.size() - 1;
	}
};