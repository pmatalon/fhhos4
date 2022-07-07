#pragma once
#include "../../Mesh/Mesh.h"
using namespace std;

template <int Dim>
class BoundaryFacePatch
{
public:
	vector<Face<Dim>*> Faces;

	// Builds the patch around a boundary face with respect to shared vertices.
	// If depth = 0, just the face f. If depth = 1, f + all its neighbours.
	// If depth = 2, + all its neighbours' neighbours. Etc.
	BoundaryFacePatch(Face<Dim>* f, vector<bool>& isFaceInAPatch, int nInteriorFaces, int depth)
	{
		set<Face<Dim>*> patch;
		set<MeshVertex<Dim>*> patchBoundaryVertices;

		AddToPatch(f, patch, patchBoundaryVertices, isFaceInAPatch, nInteriorFaces);
		
		for (int i = 0; i < depth; i++)
			AddOneLayer(patch, patchBoundaryVertices, isFaceInAPatch, nInteriorFaces);

		this->Faces = vector<Face<Dim>*>(patch.begin(), patch.end());
	}

private:
	void AddToPatch(Face<Dim>* f, set<Face<Dim>*>& patch, set<MeshVertex<Dim>*>& patchBoundaryVertices, vector<bool>& isFaceInAPatch, int nInteriorFaces)
	{
		patch.insert(f);
		isFaceInAPatch[f->Number - nInteriorFaces] = true;

		// Manage vertices
		for (Vertex* v : f->Vertices())
		{
			MeshVertex<Dim>* mv = static_cast<MeshVertex<Dim>*>(v);
			bool allFacesConnectedToThisVertexAreInPatch = true;
			for (Face<Dim>* f2 : mv->Faces)
			{
				if (f2 == f || !f2->IsDomainBoundary)
					continue;
				if (patch.find(f2) == patch.end())
				{
					allFacesConnectedToThisVertexAreInPatch = false;
					break;
				}
			}
			if (allFacesConnectedToThisVertexAreInPatch)
			{
				auto it = patchBoundaryVertices.find(mv);
				if (it != patchBoundaryVertices.end())
					patchBoundaryVertices.erase(it);
			}
			else
				patchBoundaryVertices.insert(mv);
		}
	}

	void AddOneLayer(set<Face<Dim>*>& patch, set<MeshVertex<Dim>*>& patchBoundaryVertices, vector<bool>& isFaceInAPatch, int nInteriorFaces)
	{
		vector<MeshVertex<Dim>*> currentBoundaryVertices(patchBoundaryVertices.begin(), patchBoundaryVertices.end());
		for (MeshVertex<Dim>* v : currentBoundaryVertices)
		{
			for (Face<Dim>* f : v->Faces)
			{
				if (f->IsDomainBoundary && !isFaceInAPatch[f->Number - nInteriorFaces])
					AddToPatch(f, patch, patchBoundaryVertices, isFaceInAPatch, nInteriorFaces);
			}
		}
	}

/*public:
	int Number(Face<Dim>* f) const
	{
		auto it = find(this->Faces.begin(), this->Faces.end(), f);
		if (it != this->Faces.end())
		{
			int index = it - this->Faces.begin();
			return index;
		}
		assert(false);
		return -1;
	}*/
};
