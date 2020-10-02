#pragma once
#include "../CartesianElement.h"

class ParallelepipedElement : public CartesianElement<3>
{
public:
	Face<3>* TopFace;
	Face<3>* BottomFace;
	Face<3>* FrontFace;
	Face<3>* BackFace;
	Face<3>* LeftFace;
	Face<3>* RightFace;

public:
	Vertex* BackLeftBottomCorner;
	Vertex* FrontLeftBottomCorner;
	Vertex* BackRightBottomCorner;
	Vertex* BackLeftTopCorner;
	Vertex* FrontLeftTopCorner;
	Vertex* BackRightTopCorner;
	Vertex* FrontRightBottomCorner;
	Vertex* FrontRightTopCorner;

	ParallelepipedElement(int number, Vertex* backLeftBottomCorner, Vertex* frontLeftBottomCorner, Vertex* backRightBottomCorner, Vertex* backLeftTopCorner, Vertex* frontLeftTopCorner, Vertex* backRightTopCorner, Vertex* frontRightBottomCorner, Vertex* frontRightTopCorner) :
		Element(number), 
		CartesianElement(number, backLeftBottomCorner, frontLeftBottomCorner->X - backLeftBottomCorner->X, backRightBottomCorner->Y - backLeftBottomCorner->Y, backLeftTopCorner->Z - backLeftBottomCorner->Z)
		// back  : x <-- x,    front: x <-- x + width
		// left  : y <-- y,	   right: y <-- y + width
		// bottom: z <-- z,    top  : z <-- z + width
	{
		BackLeftBottomCorner = backLeftBottomCorner;
		FrontLeftBottomCorner = frontLeftBottomCorner;
		BackRightBottomCorner = backRightBottomCorner;
		BackLeftTopCorner = backLeftTopCorner;
		FrontLeftTopCorner = frontLeftTopCorner;
		BackRightTopCorner = backRightTopCorner;
		FrontRightBottomCorner = frontRightBottomCorner;
		FrontRightTopCorner = frontRightTopCorner;
		this->SetVertices({ backLeftBottomCorner, frontLeftBottomCorner, backRightBottomCorner, backLeftTopCorner,
			frontLeftTopCorner, backRightTopCorner, frontRightBottomCorner, frontRightTopCorner});
	}

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

	DimVector<3> OuterNormalVector(Face<3>* face) const override
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