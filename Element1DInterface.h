#pragma once
#include "ElementInterface.h"

class Element1DInterface : public ElementInterface
{
public:
	double X;

	Element1DInterface(Element* element1, Element* element2) : ElementInterface(element1, element2)
	{	}

	Element1DInterface(Element* element1) : ElementInterface(element1)
	{	}
};