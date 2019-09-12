#pragma once
#include "../CartesianElement.h"

class Parallelepiped : public CartesianElement<3>
{
public:
	Face<3>* TopFace;
	Face<3>* BottomFace;
	Face<3>* FrontFace;
	Face<3>* BackFace;
	Face<3>* LeftFace;
	Face<3>* RightFace;

public:
	DomPoint BackLeftBottomCorner;
	DomPoint FrontLeftBottomCorner;
	DomPoint BackRightBottomCorner;
	DomPoint BackLeftTopCorner;

	Parallelepiped(int number, double x, double y, double z, double widthX, double widthY, double widthZ) :
		Element(number), 
		CartesianElement(number, DomPoint(x,y,z), widthX, widthY, widthZ),
		// back  : x <-- x,    front: x <-- x + width
		// left  : y <-- y,	   right: y <-- y + width
		// bottom: z <-- z,    top  : z <-- z + width
		BackLeftBottomCorner(x, y, z),
		FrontLeftBottomCorner(x + widthX, y, z),
		BackRightBottomCorner(x, y + widthY, z),
		BackLeftTopCorner(x, y, z + widthZ)
	{ }

	void SetTopFace(Face<3>* face)
	{
		this->AddFace(face);
		this->TopFace = face;
	}

	void SetBottomFace(Face<3>* face)
	{
		this->AddFace(face);
		this->BottomFace = face;
	}

	void SetFrontFace(Face<3>* face)
	{
		this->AddFace(face);
		this->FrontFace = face;
	}

	void SetBackFace(Face<3>* face)
	{
		this->AddFace(face);
		this->BackFace = face;
	}

	void SetLeftFace(Face<3>* face)
	{
		this->AddFace(face);
		this->LeftFace = face;
	}

	void SetRightFace(Face<3>* face)
	{
		this->AddFace(face);
		this->RightFace = face;
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	DimVector<3> OuterNormalVector(Face<3>* face)
	{
		DimVector<3> n;
		if (face == this->TopFace)
			n << 0, 0, 1;
		else if (face == this->BottomFace)
			n << 0, 0, -1;
		else if (face == this->FrontFace)
			n << 1, 0, 0;
		else if (face == this->BackFace)
			n << -1, 0, 0;
		else if (face == this->LeftFace)
			n << 0, -1, 0;
		else if (face == this->RightFace)
			n << 0, 1, 0;
		else
			assert(false);
		return n;
	}

};