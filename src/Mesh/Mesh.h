#pragma once
#include <vector>
#include "Element.h"
#include "Face.h"

using namespace std;

template <int Dim>
class Mesh
{
public:
	BigNumber N;
	vector<Element<Dim>*> Elements;
	vector<Face<Dim>*> Faces;
	vector<Face<Dim>*> BoundaryFaces;
	vector<Face<Dim>*> InteriorFaces;

	Mesh<Dim>* CoarserMesh = NULL;

	Mesh(BigNumber n)
	{
		this->N = n;
	}

	virtual void BuildCoarserMesh() = 0;

	virtual void Serialize(ostream& os) const
	{
		for (auto element : this->Elements)
			os << *element << endl;
		for (auto face : this->Faces)
			os << *face << endl;
	}

	friend ostream& operator<<(ostream& os, const Mesh<Dim>& s)
	{
		s.Serialize(os);
		return os;
	}

	virtual ~Mesh() 
	{
		for (size_t i = 0; i < this->Elements.size(); ++i)
			delete this->Elements[i];
		this->Elements.clear();

		for (size_t i = 0; i < this->Faces.size(); ++i)
			delete this->Faces[i];
		this->Faces.clear();

		if (CoarserMesh)
			delete CoarserMesh;
	}
};