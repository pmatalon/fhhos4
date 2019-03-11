#pragma once
#include <vector>
#include "Element.h"
#include "Face.h"

using namespace std;

template <short Dim>
class Mesh
{
public:
	BigNumber N;
	vector<Element<Dim>*> Elements;
	vector<Face<Dim>*> Faces;
	vector<Face<Dim>*> BoundaryFaces;
	vector<Face<Dim>*> InteriorFaces;

	Mesh(BigNumber n)
	{
		this->N = n;
	}

	virtual ~Mesh() 
	{
		for (size_t i = 0; i < this->Elements.size(); ++i)
			delete this->Elements[i];
		this->Elements.clear();

		for (size_t i = 0; i < this->Faces.size(); ++i)
			delete this->Faces[i];
		this->Faces.clear();
	}
};