#pragma once
#include "ElementInterface.h"

class Element1DInterface : public ElementInterface
{
public:
	double X;

	/*Element1DInterface(BigNumber number, double x, Element* element1, Element* element2) : ElementInterface(number, element1, element2)
	{	
		this->X = x;
	}

	Element1DInterface(BigNumber number, double x, Element* element1) : ElementInterface(number, element1)
	{	
		this->X = x;
	}*/

	Element1DInterface(BigNumber number, double x) : ElementInterface(number, NULL, NULL)
	{
		this->X = x;
		this->IsDomainBoundary = false;
	}
};