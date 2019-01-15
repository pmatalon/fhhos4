#pragma once
#include "Face.h"
#include "Utils.h"

class Face3D : public Face
{
public:
	bool IsInXOYPlan = false;
	bool IsInYOZPlan = false;
	bool IsInXOZPlan = false;

	Face3D(BigNumber number, Element* element1, Element* element2) : Face(number, element1, element2)
	{	}

	Face3D(BigNumber number, Element* element1) : Face(number, element1)
	{	}

	/*string ToString() override
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
	}*/
};