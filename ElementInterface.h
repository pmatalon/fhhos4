#include "Element.h"
#pragma once
class ElementInterface
{
public:
	bool IsDomainBoundary;
	Element* Element1;
	Element* Element2;

public:
	ElementInterface(Element* element1, Element* element2)
	{
		this->Element1 = element1;
		this->Element2 = element2;
	}
	ElementInterface(Element* element1)
		:ElementInterface(element1, NULL)
	{
		this->IsDomainBoundary = true;
	}

	bool IsBetween(Element* element1, Element* element2)
	{
		if (element1 == element2 && (element1 == this->Element1 || element1 == this->Element2))
			return true;
		if (element1 != this->Element1 && element1 != this->Element2)
			return false;
		if (element2 != this->Element1 && element2 != this->Element2)
			return false;
		return true;
	}

	/*double* OuterNormalVector(Element* element)
	{
		Square* square = static_cast<Square*>(element);
		return this->OuterNormalVector(square);
	}

	double* OuterNormalVector(Square* element)
	{
		if (this->)
	}*/
};