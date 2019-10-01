#pragma once
#include "../CartesianFace.h"

class RectangularFace : public CartesianFace<3>
{
public:

	RectangularFace(BigNumber number, Vertex* origin, Vertex* v2, Vertex* v3, Element<3>* element1, Element<3>* element2, CartesianShapeOrientation orientation) :
		Face(number, element1, element2), 
		CartesianFace(number, origin, firstWidth(origin, v2, orientation), secondWidth(origin, v3, orientation), element1, element2, orientation)
	{}

	RectangularFace(BigNumber number, Vertex* origin, Vertex* v2, Vertex* v3, Element<3>* element1, CartesianShapeOrientation orientation) :
		Face(number, element1), 
		CartesianFace(number, origin, firstWidth(origin, v2, orientation), secondWidth(origin, v3, orientation), element1, NULL, orientation)
	{}

private:
	static double firstWidth(Vertex* origin, Vertex* v2, CartesianShapeOrientation orientation)
	{
		if (orientation == CartesianShapeOrientation::InXOY)
			return v2->X - origin->X;
		else if (orientation == CartesianShapeOrientation::InXOZ)
			return v2->X - origin->X;
		else if (orientation == CartesianShapeOrientation::InYOZ)
			return v2->Y - origin->Y;
	}

	static double secondWidth(Vertex* origin, Vertex* v3, CartesianShapeOrientation orientation)
	{
		if (orientation == CartesianShapeOrientation::InXOY)
			return v3->Y - origin->Y;
		else if (orientation == CartesianShapeOrientation::InXOZ)
			return v3->Z - origin->Z;
		else if (orientation == CartesianShapeOrientation::InYOZ)
			return v3->Z - origin->Z;
	}
};