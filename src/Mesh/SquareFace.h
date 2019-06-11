#pragma once
#include "Face.h"
#include "CartesianFace.h"
#include "../Utils/Utils.h"
#include "Cube.h"

class SquareFace : public CartesianFace<3>
{
public:

	SquareFace(BigNumber number, DomPoint origin, double width, Element<3>* element1, Element<3>* element2, CartesianShapeOrientation orientation) :
		Face(number, element1, element2), 
		CartesianFace(number, origin, width, element1, element2, orientation)
	{ }

	SquareFace(BigNumber number, DomPoint origin, double width, Element<3>* element1, CartesianShapeOrientation orientation) :
		Face(number, element1), 
		CartesianFace(number, origin, width, element1, NULL, orientation)
	{ }
};