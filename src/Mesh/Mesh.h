#pragma once
#include <vector>
#include "Element.h"
#include "Face.h"

using namespace std;

class Mesh
{
public:
	BigNumber N;
	vector<Element*> Elements;
	vector<Face*> Faces;

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