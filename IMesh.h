#pragma once
#include <vector>
#include "Element.h"

using namespace std;

class IMesh
{
public:
	int Dim;
	BigNumber N;
	vector<Element*> Elements;
	vector<Face*> Faces;

	IMesh(int dim, BigNumber n)
	{
		this->Dim = dim;
		this->N = n;
	}

	virtual ~IMesh() 
	{
		for (size_t i = 0; i < this->Elements.size(); ++i)
			delete this->Elements[i];
		this->Elements.clear();

		for (size_t i = 0; i < this->Faces.size(); ++i)
			delete this->Faces[i];
		this->Faces.clear();
	}
};