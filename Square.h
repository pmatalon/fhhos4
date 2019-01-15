#pragma once
#include "Element.h"
#include "Face2D.h"

class Square : public Element
{
public:
	double X;
	double Y;
	double Width;

	Face2D* NorthFace;
	Face2D* SouthFace;
	Face2D* EastFace;
	Face2D* WestFace;
private:
	Square* _neighbours;
public:
	Square(int number, double x, double y, double width) : Element(number)
	{
		this->X = x;
		this->Y = y;
		this->Width = width;
	}

	StandardElementCode StdElementCode()
	{
		return StandardElementCode::Square;
	}

	void SetNorthInterface(Face2D* face)
	{
		this->Faces.push_back(face);
		this->NorthFace = face;
		face->X1 = this->X;
		face->Y1 = this->Y + this->Width;
		face->X2 = this->X + this->Width;
		face->Y2 = this->Y + this->Width;
	}

	void SetSouthInterface(Face2D* face)
	{
		this->Faces.push_back(face);
		this->SouthFace = face;
		face->X1 = this->X;
		face->Y1 = this->Y;
		face->X2 = this->X + this->Width;
		face->Y2 = this->Y;
	}

	void SetEastInterface(Face2D* face)
	{
		this->Faces.push_back(face);
		this->EastFace = face;
		face->X1 = this->X + this->Width;
		face->Y1 = this->Y;
		face->X2 = this->X + this->Width;
		face->Y2 = this->Y + this->Width;
	}

	void SetWestInterface(Face2D* face)
	{
		this->Faces.push_back(face);
		this->WestFace = face;
		face->X1 = this->X;
		face->Y1 = this->Y;
		face->X2 = this->X;
		face->Y2 = this->Y + this->Width;
	}

	double* OuterNormalVector(Face* face)
	{
		if (face == this->NorthFace)
			return new double[2]{ 0, 1 };
		if (face == this->SouthFace)
			return new double[2]{ 0, -1 };
		if (face == this->EastFace)
			return new double[2]{ 1, 0 };
		if (face == this->WestFace)
			return new double[2]{ -1, 0 };
		return NULL;
	}
};