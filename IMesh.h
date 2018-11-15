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
	vector<ElementInterface*> Interfaces;

	IMesh(int dim, BigNumber n)
	{
		this->Dim = dim;
		this->N = n;
	}
};