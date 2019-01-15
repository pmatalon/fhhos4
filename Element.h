#pragma once
#include <vector>
#include "Utils.h"
class ElementInterface;

enum class StandardElementCode
{
	None,
	Interval,
	Square,
	Cube
};

class Element 
{
public:
	BigNumber Number;
	std::vector<ElementInterface*> Interfaces;

	Element(BigNumber number)
	{
		this->Number = number;
	}

	virtual StandardElementCode StdElementCode() = 0;

	virtual double* OuterNormalVector(ElementInterface* interface) = 0;

	virtual ~Element() {}
};