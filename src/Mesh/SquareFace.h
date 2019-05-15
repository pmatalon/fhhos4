#pragma once
#include "Face.h"
#include "CartesianFace.h"
#include "../Utils/Utils.h"
#include "Cube.h"

class SquareFace : public CartesianFace<3>
{
public:

	SquareFace(BigNumber number, double width, Element<3>* element1, Element<3>* element2) : 
		Face(number, element1, element2), 
		CartesianFace(number, DomPoint(), width, element1, element2, CartesianShapeOrientation::None)
	{ }

	SquareFace(BigNumber number, double width, Element<3>* element1) : 
		Face(number, element1), 
		CartesianFace(number, DomPoint(), width, element1, NULL, CartesianShapeOrientation::None)
	{ }
};