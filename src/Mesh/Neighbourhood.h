#pragma once
#include "Element.h"
using namespace std;

template <int Dim>
class Neighbourhood
{
public:
	vector<Element<Dim>*> Elements;
	vector<Face<Dim>*> BoundaryFaces;
	vector<Face<Dim>*> InteriorFaces;

	Neighbourhood(Element<Dim>* e)
	{
		set<Element<Dim>*> neighbours;
		neighbours.insert(e);
		for (Face<Dim>* f : e->Faces)
		{
			if (f->IsDomainBoundary)
			{
				BoundaryFaces.push_back(f);
				continue;
			}
			neighbours.insert(f->GetNeighbour(e));
			InteriorFaces.push_back(f);
		}

		Elements = vector<Element<Dim>*>(neighbours.begin(), neighbours.end());

		for (Element<Dim>* n : Elements)
		{
			if (n == e)
				continue;
			for (Face<Dim>* f : n->Faces)
			{
				if (!f->IsIn(InteriorFaces))
					BoundaryFaces.push_back(f);
			}
		}
	}

	int ElementNumber(Element<Dim>* e) const
	{
		auto it = find(Elements.begin(), Elements.end(), e);
		if (it != Elements.end())
		{
			int index = it - Elements.begin();
			return index;
		}
		assert(false);
		return -1;
	}

	int FaceNumber(Face<Dim>* f) const
	{
		auto it = find(InteriorFaces.begin(), InteriorFaces.end(), f);
		if (it != InteriorFaces.end())
		{
			int index = it - InteriorFaces.begin();
			return index;
		}
		int boundaryNumber = BoundaryFaceNumber(f);
		if (boundaryNumber != -1)
			return InteriorFaces.size() + boundaryNumber;
		assert(false);
		return -1;
	}

	int BoundaryFaceNumber(Face<Dim>* f) const
	{
		auto it = find(BoundaryFaces.begin(), BoundaryFaces.end(), f);
		if (it != BoundaryFaces.end())
		{
			int index = it - BoundaryFaces.begin();
			return index;
		}
		assert(false);
		return -1;
	}
};
