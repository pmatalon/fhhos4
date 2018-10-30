#pragma once
#include "ElementInterface.h"
#include "Utils.h"
#include <functional>
class Element2DInterface : public ElementInterface
{
public:
	double X1;
	double Y1;
	double X2;
	double Y2;

	Element2DInterface(Element* element1, Element* element2) : ElementInterface(element1, element2)
	{	}

	Element2DInterface(Element* element1) : ElementInterface(element1)
	{	}

	double Integrate(std::function<double(double, double)> func)
	{
		return Utils::Integral(func, this->X1, this->X2, this->Y1, this->Y2);
	}
};