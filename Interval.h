#pragma once
#include "Element.h"
#include "Element1DInterface.h"

class Interval : public Element
{
public:
	double A;
	double B;

	Element1DInterface* LeftInterface;
	Element1DInterface* RightInterface;

	Interval(BigNumber number, double a, double b) : Element(number)
	{
		this->A = a;
		this->B = b;
	}

	void SetLeftInterface(Element1DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->LeftInterface = interface;
		interface->X = this->A;
	}

	void SetRightInterface(Element1DInterface* interface)
	{
		this->Interfaces.push_back(interface);
		this->RightInterface = interface;
		interface->X = this->B;
	}

	double* OuterNormalVector(ElementInterface* interface)
	{
		if (interface == this->LeftInterface)
			return new double[1]{ -1 };
		else if(interface == this->RightInterface)
			return new double[1]{ 1 };
		return NULL;
	}
};