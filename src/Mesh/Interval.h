#pragma once
#include "CartesianElement.h"

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

	vector<double> OuterNormalVector(Face<1>* interface)
	{
		if (interface == this->Left)
			return vector<double>{ -1 };
		else if(interface == this->Right)
			return vector<double>{ 1 };
		assert(false);
	}
};