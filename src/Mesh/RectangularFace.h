#pragma once
#include "Face.h"
#include "CartesianFace.h"
#include "../Utils/Utils.h"

class RectangularFace : public CartesianFace<3>
{
public:

	RectangularFace(BigNumber number, DomPoint origin, double firstWidth, double secondWidth, Element<3>* element1, Element<3>* element2, CartesianShapeOrientation orientation) :
		Face(number, element1, element2), 
		CartesianFace(number, origin, firstWidth, secondWidth, element1, element2, orientation)
	{ }

	RectangularFace(BigNumber number, DomPoint origin, double firstWidth, double secondWidth, Element<3>* element1, CartesianShapeOrientation orientation) :
		Face(number, element1), 
		CartesianFace(number, origin, firstWidth, secondWidth, element1, NULL, orientation)
	{ }
};