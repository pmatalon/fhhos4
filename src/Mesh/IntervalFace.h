#pragma once
#include "Face.h"
#include "CartesianFace.h"
#include "../Utils/Utils.h"
#include "../DG/Poisson_DG_Face.h"

class IntervalFace : public CartesianFace<2>
{
public:

	IntervalFace(BigNumber number, DomPoint origin, double length, Element<2>* element1, Element<2>* element2, CartesianShapeOrientation orientation) : 
		Face(number, element1, element2), 
		CartesianFace(number, origin, length, element1, element2, orientation)
	{ }

	IntervalFace(BigNumber number, DomPoint origin, double length, Element<2>* element1, CartesianShapeOrientation orientation) :
		Face(number, element1), 
		CartesianFace(number, origin, length, element1, NULL, orientation)
	{ }

	friend ostream& operator<<(ostream& os, const IntervalFace& face)
	{
		os << face.Number << "\t(" << face.CartesianShape::Origin.X << ", " << face.CartesianShape::Origin.Y << ")\t" << (face.CartesianShape::Orientation == CartesianShapeOrientation::Horizontal ? "horizontal" : "vertical");
		return os;
	}
};