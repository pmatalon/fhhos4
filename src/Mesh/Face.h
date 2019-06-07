#include "Element.h"
#pragma once

template <int Dim>
class Face
{
public:
	BigNumber Number;
	bool IsDomainBoundary;
	Element<Dim>* Element1;
	Element<Dim>* Element2;

public:
	Face() { assert(false); }
	Face(BigNumber number, Element<Dim>* element1, Element<Dim>* element2)
	{
		this->Number = number;
		this->Element1 = element1;
		this->Element2 = element2;
		this->IsDomainBoundary = element2 == NULL;
	}
	Face(BigNumber number, Element<Dim>* element1)
		:Face(number, element1, NULL)
	{
		this->IsDomainBoundary = true;
	}

	bool IsBetween(Element<Dim>* element1, Element<Dim>* element2)
	{
		if (element1 == element2 && (element1 == this->Element1 || element1 == this->Element2))
			return true;
		if (element1 != this->Element1 && element1 != this->Element2)
			return false;
		if (element2 != this->Element1 && element2 != this->Element2)
			return false;
		return true;
	}

	Element<Dim>* GetNeighbour(Element<Dim>* element)
	{
		if (element == this->Element1)
			return this->Element2;
		else if (element == this->Element2)
			return this->Element1;
		return NULL;
	}

	virtual double GetDiameter() = 0;
	virtual double Measure() = 0;
	virtual DomPoint ConvertToDomain(RefPoint refPoint) = 0;
	virtual double ComputeIntegral(function<double(RefPoint)> func) = 0;
	virtual double ComputeIntegral(function<double(RefPoint)> func, int polynomialDegree) = 0;
	
	virtual string ToString()
	{
		return "Interface " + to_string(this->Number);
	}

	virtual ~Face() {}
};