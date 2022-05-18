#pragma once
#include "../BasisFunction.h"
#include <math.h>
using namespace std;

class Monomial1D : public IBasisFunction1D
{
public:
	int Degree;

	Monomial1D(int degree)
	{
		this->LocalNumber = degree;
		this->Degree = degree;
	}

	int GetDegree() const
	{
		return this->Degree;
	}

	double Eval(double x)
	{
		this->TestIsInReferenceInterval(x);
		return pow(x, this->Degree);
	}

	double EvalDerivative(double x)
	{
		this->TestIsInReferenceInterval(x);
		if (this->Degree == 0)
			return 0;
		return this->Degree*pow(x, this->Degree - 1);
	}

	string ToString()
	{
		return this->ToString("X");
	}

	string ToString(string var)
	{
		if (this->Degree == 0)
			return "1";
		if (this->Degree == 1)
			return var;
		return var + "^" + std::to_string(this->Degree);
	}
};