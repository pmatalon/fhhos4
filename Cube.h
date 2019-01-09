#pragma once
#include "Element.h"
#include "Element3DInterface.h"

class Cube : public Element
{
public:
	double X;
	double Y;
	double Z;
	double Width;

	Element3DInterface* TopInterface;
	Element3DInterface* BottomInterface;
	Element3DInterface* FrontInterface;
	Element3DInterface* BackInterface;
	Element3DInterface* LeftInterface;
	Element3DInterface* RightInterface;
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

	void SetTopInterface(Element3DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->TopInterface = interface;
		interface->IsInXOYPlan = true;
	}

	void SetBottomInterface(Element3DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->BottomInterface = interface;
		interface->IsInXOYPlan = true;
	}

	void SetFrontInterface(Element3DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->FrontInterface = interface;
		interface->IsInXOZPlan = true;
	}

	void SetBackInterface(Element3DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->BackInterface = interface;
		interface->IsInXOZPlan = true;
	}

	void SetLeftInterface(Element3DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->LeftInterface = interface;
		interface->IsInYOZPlan = true;
	}

	void SetRightInterface(Element3DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->RightInterface = interface;
		interface->IsInYOZPlan = true;
	}

	double* OuterNormalVector(ElementInterface* interface)
	{
		if (interface == this->TopInterface)
			return new double[3]{ 0, 0, 1 };
		if (interface == this->BottomInterface)
			return new double[3]{ 0, 0, -1 };
		if (interface == this->FrontInterface)
			return new double[3]{ 0, -1, 0 };
		if (interface == this->BackInterface)
			return new double[3]{ 0, 1, 0 };
		if (interface == this->LeftInterface)
			return new double[3]{ -1, 0, 0 };
		if (interface == this->RightInterface)
			return new double[3]{ 1, 0, 0 };
		return NULL;
	}
};