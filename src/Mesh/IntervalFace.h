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
};