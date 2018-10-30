#pragma once
#include "BasisFunction1D.h"
#include <math.h>

class ReverseMonomial1D : public BasisFunction1D
{
public:
	int Degree;

	ReverseMonomial1D(int degree)
	{
		this->Degree = degree;
	}

	double Eval(double x)
	{
		return pow(1 - x, this->Degree);
	}

	double EvalGrad(double x)
	{
		if (this->Degree == 0)
			return 0;
		return -this->Degree*pow(1 - x, this->Degree - 1);
	}
};
