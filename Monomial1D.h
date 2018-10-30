#pragma once
#include "BasisFunction1D.h"
#include <math.h>

class Monomial1D : public BasisFunction1D
{
public:
	int Degree;

	Monomial1D(int degree)
	{
		this->Degree = degree;
	}

	double Eval(double x)
	{
		return pow(x, this->Degree);
	}

	double EvalGrad(double x)
	{
		if (this->Degree == 0)
			return 0;
		return this->Degree*pow(x, this->Degree - 1);
	}
};