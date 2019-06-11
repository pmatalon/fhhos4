#pragma once
#include "CartesianElement.h"
#include "Edge.h"
using namespace std;

class Rectangle : public CartesianElement<2>
{
public:
	Face<2>* NorthFace;
	Face<2>* SouthFace;
	Face<2>* EastFace;
	Face<2>* WestFace;

	DomPoint BottomLeftCorner;
	DomPoint TopLeftCorner;
	DomPoint TopRightCorner;
	DomPoint BottomRightCorner;

	Rectangle(int number, double x, double y, double widthX, double widthY) :
		Element(number), 
		CartesianElement(number, DomPoint(x, y), widthX, widthY), 
		BottomLeftCorner(x, y),
		TopLeftCorner(x, y + widthY),
		TopRightCorner(x + widthX, y + widthY),
		BottomRightCorner(x + widthX, y)
	{}

	void SetNorthInterface(Edge* face)
	{
		this->AddFace(face);
		this->NorthFace = face;
	}

	void SetSouthInterface(Edge* face)
	{
		this->AddFace(face);
		this->SouthFace = face;
	}

	void SetEastInterface(Edge* face)
	{
		this->AddFace(face);
		this->EastFace = face;
	}

	void SetWestInterface(Edge* face)
	{
		this->AddFace(face);
		this->WestFace = face;
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	vector<double> OuterNormalVector(Face<2>* face)
	{
		if (face == this->NorthFace)
			return vector<double>{ 0, 1 };
		if (face == this->SouthFace)
			return vector<double>{ 0, -1 };
		if (face == this->EastFace)
			return vector<double>{ 1, 0 };
		if (face == this->WestFace)
			return vector<double>{ -1, 0 };
		assert(false);
	}

};