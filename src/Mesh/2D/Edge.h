#pragma once
#include "../CartesianFace.h"

class Edge : public CartesianFace<2>
{
public:

	Edge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1, Element<2>* element2, CartesianShapeOrientation orientation) : 
		Face(number, element1, element2), 
		CartesianFace(number, v1, orientation == CartesianShapeOrientation::Horizontal ? v2->X - v1->X : v2->Y - v1->Y, element1, element2, orientation)
	{ 
	}

	Edge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1, CartesianShapeOrientation orientation) :
		Face(number, element1), 
		CartesianFace(number, v1, orientation == CartesianShapeOrientation::Horizontal ? v2->X - v1->X : v2->Y - v1->Y, element1, NULL, orientation)
	{
	}
};