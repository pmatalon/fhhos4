#pragma once
#include "../CartesianFace.h"

class Edge : public CartesianFace<2>
{
public:

	Edge(BigNumber number, DomPoint origin, double length, Element<2>* element1, Element<2>* element2, CartesianShapeOrientation orientation) : 
		Face(number, element1, element2), 
		CartesianFace(number, origin, length, element1, element2, orientation)
	{ }

	Edge(BigNumber number, DomPoint origin, double length, Element<2>* element1, CartesianShapeOrientation orientation) :
		Face(number, element1), 
		CartesianFace(number, origin, length, element1, NULL, orientation)
	{ }
};