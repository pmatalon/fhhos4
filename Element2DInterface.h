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

	Element2DInterface(BigNumber number, Element* element1, Element* element2) : ElementInterface(number, element1, element2)
	{	}

	Element2DInterface(BigNumber number, Element* element1) : ElementInterface(number, element1)
	{	}

	/*double Integrate(std::function<double(double, double)> func)
	{
		return Utils::Integral(func, this->X1, this->X2, this->Y1, this->Y2);
	}*/

	bool IsVertical() {	return this->X1 == this->X2; }
	bool IsHorizontal()	{ return this->Y1 == this->Y2; }

	string ToString() override
	{
		string s = "Interface " + std::to_string(this->Number);
		if (this->IsVertical())
			s += " (vertical)";
		else if (this->IsHorizontal())
			s += " (horizontal)";
		if (this->IsDomainBoundary)
			s += " on element " + std::to_string(this->Element1->Number) + " (boundary)";
		else
			s += " between element " + std::to_string(this->Element1->Number) + " and element " + std::to_string(this->Element2->Number);
		return s;
	}
};