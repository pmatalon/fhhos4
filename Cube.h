#pragma once
#include "Element.h"
#include "Face3D.h"

class Cube : public Element
{
public:
	double X;
	double Y;
	double Z;
	double Width;

	Face3D* TopFace;
	Face3D* BottomFace;
	Face3D* FrontFace;
	Face3D* BackFace;
	Face3D* LeftFace;
	Face3D* RightFace;
private:
	Cube* _neighbours;
public:
	Cube(int number, double x, double y, double z, double width) : Element(number)
	{
		this->X = x;
		this->Y = y;
		this->Z = z;
		this->Width = width;
	}

	StandardElementCode StdElementCode()
	{
		return StandardElementCode::Cube;
	}

	void SetTopInterface(Face3D* interface)
	{
		this->Faces.push_back(interface);
		this->TopFace = interface;
		interface->IsInXOYPlan = true;
	}

	void SetBottomInterface(Face3D* interface)
	{
		this->Faces.push_back(interface);
		this->BottomFace = interface;
		interface->IsInXOYPlan = true;
	}

	void SetFrontInterface(Face3D* interface)
	{
		this->Faces.push_back(interface);
		this->FrontFace = interface;
		interface->IsInXOZPlan = true;
	}

	void SetBackInterface(Face3D* interface)
	{
		this->Faces.push_back(interface);
		this->BackFace = interface;
		interface->IsInXOZPlan = true;
	}

	void SetLeftInterface(Face3D* interface)
	{
		this->Faces.push_back(interface);
		this->LeftFace = interface;
		interface->IsInYOZPlan = true;
	}

	void SetRightInterface(Face3D* interface)
	{
		this->Faces.push_back(interface);
		this->RightFace = interface;
		interface->IsInYOZPlan = true;
	}

	double* OuterNormalVector(Face* interface)
	{
		if (interface == this->TopFace)
			return new double[3]{ 0, 0, 1 };
		if (interface == this->BottomFace)
			return new double[3]{ 0, 0, -1 };
		if (interface == this->FrontFace)
			return new double[3]{ 0, -1, 0 };
		if (interface == this->BackFace)
			return new double[3]{ 0, 1, 0 };
		if (interface == this->LeftFace)
			return new double[3]{ -1, 0, 0 };
		if (interface == this->RightFace)
			return new double[3]{ 1, 0, 0 };
		return NULL;
	}
};