#pragma once
#include "../CartesianFace.h"

class CartesianEdge : public CartesianFace<2>
{
public:
	CartesianEdge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1, Element<2>* element2, CartesianShapeOrientation orientation) : 
		Face(number, element1, element2), 
		CartesianFace(number, v1, orientation == CartesianShapeOrientation::Horizontal ? v2->X - v1->X : v2->Y - v1->Y, element1, element2, orientation)
	{
		this->_shape.SetVertices(vector<Vertex*> {v1, v2});
	}

	CartesianEdge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1, CartesianShapeOrientation orientation) :
		Face(number, element1), 
		CartesianFace(number, v1, orientation == CartesianShapeOrientation::Horizontal ? v2->X - v1->X : v2->Y - v1->Y, element1, NULL, orientation)
	{
		this->_shape.SetVertices(vector<Vertex*> {v1, v2});
	}
};