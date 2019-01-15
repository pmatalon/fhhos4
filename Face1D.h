#pragma once
#include "Face.h"

class Face1D : public Face
{
public:
	double X;

	/*Face1D(BigNumber number, double x, Element* element1, Element* element2) : Face(number, element1, element2)
	{	
		this->X = x;
	}

	Face1D(BigNumber number, double x, Element* element1) : Face(number, element1)
	{	
		this->X = x;
	}*/

	Face1D(BigNumber number, double x) : Face(number, NULL, NULL)
	{
		this->X = x;
		this->IsDomainBoundary = false;
	}
};