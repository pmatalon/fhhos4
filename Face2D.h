#pragma once
#include "Face.h"
#include "Utils.h"

class Face2D : public Face
{
public:
	double X1;
	double Y1;
	double X2;
	double Y2;

	Face2D(BigNumber number, Element* element1, Element* element2) : Face(number, element1, element2)
	{	}

	Face2D(BigNumber number, Element* element1) : Face(number, element1)
	{	}

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