#pragma once
#include "IBasisFunction1D.h"
#include <math.h>
using namespace std;

class Monomial1D : public IBasisFunction1D
{
public:
	int Degree;

	static string Code() { return "monomials"; };

	Monomial1D(int degree)
	{
		this->Degree = degree;
	}

	RefInterval ReferenceInterval() { return RefInterval::MinusOne_One(); }

	int GetDegree()
	{
		return this->Degree;
	}

	double Eval(double x)
	{
		return pow(x, this->Degree);
	}

	double EvalDerivative(double x)
	{
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

class MonomialBasis1D : public FunctionalBasisWithNumbers
{
private:
	int _maxPolynomialDegree;

public:
	MonomialBasis1D(int maxPolynomialDegree)
		:FunctionalBasisWithNumbers()
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		for (int i = 0; i <= maxPolynomialDegree; i++)
			this->_localFunctions[i] = new Monomial1D(i);
	}

	int GetDegree()
	{
		return this->_maxPolynomialDegree;
	}

	string Name()
	{
		return Monomial1D::Code() + "_p" + std::to_string(this->_maxPolynomialDegree);
	}
};