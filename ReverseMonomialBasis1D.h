/*#pragma once
#include "FunctionalBasis1D.h"
#include "IBasisFunction1D.h"
#include "IPolynomialFunction.h"

using namespace std;

class ReverseMonomial1D : public IBasisFunction1D, public IPolynomialFunction
{
public:
	int Degree;

	ReverseMonomial1D(int degree)
	{
		this->Degree = degree;
	}

	int GetDegree()
	{
		return this->Degree;
	}

	double Eval(double x)
	{
		return pow(1 - x, this->Degree);
	}

	double EvalDerivative(double x)
	{
		if (this->Degree == 0)
			return 0;
		return -this->Degree*pow(1 - x, this->Degree - 1);
	}
};

class ReverseMonomialBasis1D : public FunctionalBasis1D
{
private:
	int _maxPolynomialDegree;

public:
	ReverseMonomialBasis1D(int maxPolynomialDegree, CartesianGrid1D* grid, int penalizationCoefficient, function<double(double)> sourceFunction)
		:FunctionalBasis1D(grid, penalizationCoefficient, sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		for (int i = 0; i <= maxPolynomialDegree; i++)
			this->_localFunctions[i] = new ReverseMonomial1D(i);
	}

	string Name()
	{
		return "reversemonomials_p" + std::to_string(this->_maxPolynomialDegree);
	}
};*/