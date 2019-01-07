#pragma once
#include "Element.h"
#include "Element2DInterface.h"

class Square : public Element
{
public:
	double X;
	double Y;
	double Width;

	Element2DInterface* NorthInterface;
	Element2DInterface* SouthInterface;
	Element2DInterface* EastInterface;
	Element2DInterface* WestInterface;
private:
	Square* _neighbours;
public:
	Square(int number, double x, double y, double width) : Element(number)
	{
		this->X = x;
		this->Y = y;
		this->Width = width;
	}

	void SetNorthInterface(Element2DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->NorthInterface = interface;
		interface->X1 = this->X;
		interface->Y1 = this->Y + this->Width;
		interface->X2 = this->X + this->Width;
		interface->Y2 = this->Y + this->Width;
	}

	void SetSouthInterface(Element2DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->SouthInterface = interface;
		interface->X1 = this->X;
		interface->Y1 = this->Y;
		interface->X2 = this->X + this->Width;
		interface->Y2 = this->Y;
	}

	void SetEastInterface(Element2DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->EastInterface = interface;
		interface->X1 = this->X + this->Width;
		interface->Y1 = this->Y;
		interface->X2 = this->X + this->Width;
		interface->Y2 = this->Y + this->Width;
	}

	void SetWestInterface(Element2DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->WestInterface = interface;
		interface->X1 = this->X;
		interface->Y1 = this->Y;
		interface->X2 = this->X;
		interface->Y2 = this->Y + this->Width;
	}

	double* OuterNormalVector(ElementInterface* interface)
	{
		if (interface == this->NorthInterface)
			return new double[2]{ 0, 1 };
		if (interface == this->SouthInterface)
			return new double[2]{ 0, -1 };
		if (interface == this->EastInterface)
			return new double[2]{ 1, 0 };
		if (interface == this->WestInterface)
			return new double[2]{ -1, 0 };
		return NULL;
	}
};