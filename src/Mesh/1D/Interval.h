#pragma once
#include "../CartesianElement.h"

class Interval : public CartesianElement<1>
{
public:
	Face<1>* Left;
	Face<1>* Right;

	Interval(BigNumber number, Vertex* left, Vertex* right) :
		Element(number),
		CartesianElement(number, left, right->X - left->X)
	{ }

	inline double Width()
	{
		return this->_shape.WidthX;
	}
	
	void SetLeftInterface(Face<1>* face)
	{
		this->AddFace(face);
		this->Left = face;
	}

	void SetRightInterface(Face<1>* face)
	{
		this->AddFace(face);
		this->Right = face;
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