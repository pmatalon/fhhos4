#pragma once
#include "../CartesianElement.h"

class Interval : public CartesianElement<1>
{
public:
	Face<1>* Left;
	Face<1>* Right;

	Interval(BigNumber number, double a, double b, Face<1>* left, Face<1>* right) : 
		Element(number),
		CartesianElement(number, DomPoint(a), b-a)
	{
		this->AddFace(left);
		this->AddFace(right);
		this->Left = left;
		this->Right = right;
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	DimVector<1> OuterNormalVector(Face<1>* interface)
	{
		DimVector<1> n;
		if (interface == this->Left)
			n << -1;
		else if (interface == this->Right)
			n << 1;
		else
			assert(false);
		return n;
	}
};