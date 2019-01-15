#pragma once
#include "Element.h"
#include "Face1D.h"

class Interval : public Element
{
public:
	double A;
	double B;

	Face1D* Left;
	Face1D* Right;

	Interval(BigNumber number, Face1D* left, Face1D* right) : Element(number)
	{
		this->A = left->X;
		this->B = right->X;
		this->Faces.push_back(left);
		this->Faces.push_back(right);
		this->Left = left;
		this->Right = right;
	}

	StandardElementCode StdElementCode()
	{
		return StandardElementCode::Interval;
	}

	double* OuterNormalVector(Face* interface)
	{
		if (interface == this->Left)
			return new double[1]{ -1 };
		else if(interface == this->Right)
			return new double[1]{ 1 };
		return NULL;
	}
};