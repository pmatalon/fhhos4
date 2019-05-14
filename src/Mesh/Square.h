#pragma once
#include "CartesianElement.h"
#include "IntervalFace.h"
#include "../Utils/SourceFunction.h"
#include <assert.h>
using namespace std;

class Square : public CartesianElement<2>
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

	Square(int number, double x, double y, double width) : 
		Element(number), 
		CartesianElement(number, DomPoint(x, y), width), 
		BottomLeftCorner(x, y),
		TopLeftCorner(x, y + width),
		TopRightCorner(x + width, y + width),
		BottomRightCorner(x + width, y)
	{
	}

	void SetNorthInterface(IntervalFace* face)
	{
		this->AddFace(face);
		this->NorthFace = face;
	}

	void SetSouthInterface(IntervalFace* face)
	{
		this->AddFace(face);
		this->SouthFace = face;
	}

	void SetEastInterface(IntervalFace* face)
	{
		this->AddFace(face);
		this->EastFace = face;
	}

	void SetWestInterface(IntervalFace* face)
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