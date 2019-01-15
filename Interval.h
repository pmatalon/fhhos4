#pragma once
#include "Element.h"
#include "Element1DInterface.h"

class Interval : public Element
{
public:
	double A;
	double B;

	Element1DInterface* Left;
	Element1DInterface* Right;

	Interval(BigNumber number, Element1DInterface* left, Element1DInterface* right) : Element(number)
	{
		this->A = left->X;
		this->B = right->X;
		this->Interfaces.push_back(left);
		this->Interfaces.push_back(right);
		this->Left = left;
		this->Right = right;
	}

	StandardElementCode StdElementCode()
	{
		return StandardElementCode::Interval;
	}

	double* OuterNormalVector(ElementInterface* interface)
	{
		if (interface == this->Left)
			return new double[1]{ -1 };
		else if(interface == this->Right)
			return new double[1]{ 1 };
		return NULL;
	}
};