#pragma once
#include "IBasisFunction.h"
#include <string>
#include <math.h>
using namespace std;

class Monomial1D : public IBasisFunction1D
{
public:
	int Degree;

	static string Code() { return "monomials"; };

	Monomial1D(int degree)
	{
		this->LocalNumber = degree;
		this->Degree = degree;
	}

	DefInterval DefinitionInterval() { return DefInterval::MinusOne_One(); }

	int GetDegree()
	{
		return this->Degree;
	}

	double Eval(double x)
	{
		assert(x >= -1 && x <= 1);
		return pow(x, this->Degree);
	}

	double EvalDerivative(double x)
	{
		assert(x >= -1 && x <= 1);
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