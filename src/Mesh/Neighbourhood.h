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

	// Builds the neighbourhood of an element with respect to shared faces.
	// If depth = 0, just the element e. If depth = 1, e + all its neighbours.
	// If depth = 2, + all its neighbours' neighbours. Etc.
	Neighbourhood(Element<Dim>* e, int depth = 1)
	{
		set<Element<Dim>*> patch;
		set<Face<Dim>*> boundaryFaces;
		set<Face<Dim>*> interiorFaces;

		patch.insert(e);

		AddNeighboursOf(e, depth, patch, boundaryFaces, interiorFaces);

		Elements = vector<Element<Dim>*>(patch.begin(), patch.end());
		BoundaryFaces = vector<Face<Dim>*>(boundaryFaces.begin(), boundaryFaces.end());
		InteriorFaces = vector<Face<Dim>*>(interiorFaces.begin(), interiorFaces.end());
	}

private:
	void AddNeighboursOf(Element<Dim>* e, int depth, set<Element<Dim>*>& patch, set<Face<Dim>*>& boundaryFaces, set<Face<Dim>*>& interiorFaces)
	{
		if (depth == 0)
		{
			// Add no neighbour.
			// All the faces that are not interior to the patch is a boundary of the patch
			for (Face<Dim>* f : e->Faces)
			{
				if (interiorFaces.find(f) == interiorFaces.end())
					boundaryFaces.insert(f);
			}
			return;
		}

		for (Face<Dim>* f : e->Faces)
		{
			if (f->IsDomainBoundary)
			{
				boundaryFaces.insert(f);
				continue;
			}
			
			Element<Dim>* n = f->GetNeighbour(e);
			
			interiorFaces.insert(f);
			auto it = boundaryFaces.find(f);
			if (it != boundaryFaces.end())
				boundaryFaces.erase(it);

			// If neighbour not already included
			if (patch.find(n) == patch.end())
			{
				// Add neighbour
				patch.insert(n);
				// Add its own neighbours
				AddNeighboursOf(n, depth - 1, patch, boundaryFaces, interiorFaces);
			}
		}
	}

public:
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
