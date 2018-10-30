#pragma once
#include "BasisFunction2D.h"
#include <math.h>

class Monomial2D : public BasisFunction2D
{
public:
	int DegreeX;
	int DegreeY;

	Monomial2D(int degreeX, int degreeY)
	{
		this->DegreeX = degreeX;
		this->DegreeY = degreeY;
	}

	double Eval(double x, double y)
	{
		return pow(x, this->DegreeX)*pow(y, this->DegreeY);
	}

	double EvalGradX(double x, double y)
	{
		if (this->DegreeX == 0)
			return 0;
		return this->DegreeX*pow(x, this->DegreeX - 1)*pow(y, this->DegreeY);
	}

	double EvalGradY(double x, double y)
	{
		if (this->DegreeY == 0)
			return 0;
		return this->DegreeY*pow(y, this->DegreeY - 1)*pow(x, this->DegreeX);
	}
};