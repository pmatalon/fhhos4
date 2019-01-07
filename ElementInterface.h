#include "Element.h"
#pragma once
class ElementInterface
{
public:
	BigNumber Number;
	bool IsDomainBoundary;
	Element* Element1;
	Element* Element2;

public:
	ElementInterface(BigNumber number, Element* element1, Element* element2)
	{
		this->Number = number;
		this->Element1 = element1;
		this->Element2 = element2;
		this->IsDomainBoundary = false;
	}
	ElementInterface(BigNumber number, Element* element1)
		:ElementInterface(number, element1, NULL)
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

	virtual string ToString()
	{
		return "Interface " + this->Number;
	}

	virtual ~ElementInterface() {}
};