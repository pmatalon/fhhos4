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

	// Builds the neighbourhood of an element with respect to shared vertices.
	// If depth = 0, just the element e. If depth = 1, e + all its neighbours.
	// If depth = 2, + all its neighbours' neighbours. Etc.
	Neighbourhood(Element<Dim>* e, int depth = 1)
	{
		set<Element<Dim>*> patch;
		set<MeshVertex<Dim>*> boundaryVertices;
		set<Face<Dim>*> boundaryFaces;
		set<Face<Dim>*> interiorFaces;

		AddToPatch(e, patch, boundaryVertices, boundaryFaces, interiorFaces);
		
		for (int i = 0; i < depth; i++)
			AddOneLayer(patch, boundaryVertices, boundaryFaces, interiorFaces);

		Elements = vector<Element<Dim>*>(patch.begin(), patch.end());
		BoundaryFaces = vector<Face<Dim>*>(boundaryFaces.begin(), boundaryFaces.end());
		InteriorFaces = vector<Face<Dim>*>(interiorFaces.begin(), interiorFaces.end());
	}

private:
	void AddToPatch(Element<Dim>* e, set<Element<Dim>*>& patch, set<MeshVertex<Dim>*>& boundaryVertices, set<Face<Dim>*>& boundaryFaces, set<Face<Dim>*>& interiorFaces)
	{
		patch.insert(e);

		// Manage faces
		for (Face<Dim>* f : e->Faces)
		{
			if (f->IsDomainBoundary || patch.find(f->GetNeighbour(e)) == patch.end())
				boundaryFaces.insert(f);
			else
			{
				interiorFaces.insert(f);
				auto it = boundaryFaces.find(f);
				if (it != boundaryFaces.end())
					boundaryFaces.erase(it);
			}
		}

		// Manage vertices
		for (Vertex* v : e->Vertices())
		{
			MeshVertex<Dim>* mv = static_cast<MeshVertex<Dim>*>(v);
			bool allInPatch = true;
			for (Element<Dim>* ve : mv->Elements)
			{
				if (ve == e)
					continue;
				if (patch.find(ve) == patch.end())
				{
					allInPatch = false;
					break;
				}
			}
			if (allInPatch)
			{
				auto it = boundaryVertices.find(mv);
				if (it != boundaryVertices.end())
					boundaryVertices.erase(it);
			}
			else
				boundaryVertices.insert(mv);
		}
	}

	void AddOneLayer(set<Element<Dim>*>& patch, set<MeshVertex<Dim>*>& boundaryVertices, set<Face<Dim>*>& boundaryFaces, set<Face<Dim>*>& interiorFaces)
	{
		vector<MeshVertex<Dim>*> currentBoundaryVertices(boundaryVertices.begin(), boundaryVertices.end());
		for (MeshVertex<Dim>* v : currentBoundaryVertices)
		{
			for (Element<Dim>* e : v->Elements)
			{
				if (patch.find(e) == patch.end())
					AddToPatch(e, patch, boundaryVertices, boundaryFaces, interiorFaces);
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
