#pragma once
#include "../CartesianElement.h"
#include "CartesianEdge.h"
using namespace std;

class RectangularElement : public CartesianElement<2>
{
public:
	Face<2>* NorthFace;
	Face<2>* SouthFace;
	Face<2>* EastFace;
	Face<2>* WestFace;

	Vertex* BottomLeftCorner;
	Vertex* TopLeftCorner;
	Vertex* TopRightCorner;
	Vertex* BottomRightCorner;

	RectangularElement(int number, Vertex* bottomLeftCorner, Vertex* topLeftCorner, Vertex* topRightCorner, Vertex* bottomRightCorner) :
		Element(number), 
		CartesianElement(number, bottomLeftCorner, bottomRightCorner->X - bottomLeftCorner->X, topLeftCorner->Y - bottomLeftCorner->Y)
	{
		BottomLeftCorner = bottomLeftCorner;
		TopLeftCorner = topLeftCorner;
		TopRightCorner = topRightCorner;
		BottomRightCorner = bottomRightCorner;
		this->_shape->SetVertices(vector<Vertex*> {bottomLeftCorner, bottomRightCorner, topRightCorner, topLeftCorner});
	}

	void SetNorthInterface(CartesianEdge* face)
	{
		this->AddFace(face);
		this->NorthFace = face;
	}

	void SetSouthInterface(CartesianEdge* face)
	{
		this->AddFace(face);
		this->SouthFace = face;
	}

	void SetEastInterface(CartesianEdge* face)
	{
		this->AddFace(face);
		this->EastFace = face;
	}

	void SetWestInterface(CartesianEdge* face)
	{
		this->AddFace(face);
		this->WestFace = face;
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	DimVector<2> OuterNormalVector(Face<2>* face) const
	{
		DimVector<2> n;
		if (face == this->NorthFace)
			n << 0, 1;
		else if (face == this->SouthFace)
			n << 0, -1;
		else if (face == this->EastFace)
			n << 1, 0;
		else if (face == this->WestFace)
			n << -1, 0;
		else
			assert(false);
		return n;
	}

};