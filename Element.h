//#include "ElementInterface.h"
#include <vector>

class ElementInterface;

typedef unsigned int BigNumber;

#pragma once
class Element 
{
public:
	BigNumber Number;
	std::vector<ElementInterface*> Interfaces;

	Element(BigNumber number)
	{
		this->Number = number;
	}

	virtual double* OuterNormalVector(ElementInterface* interface) = 0;

	virtual ~Element() {}
};