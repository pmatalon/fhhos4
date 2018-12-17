#pragma once
#include "FunctionalBasis1D.h"
#include "FunctionalGlobalBasis1D.h"
#include "IBasisFunction1D.h"
#include <math.h>
using namespace std;

class Monomial1D : public IBasisFunction1D
{
public:
	int Degree;

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
		if (this->Degree == 0)
			return "1";
		if (this->Degree == 1)
			return "X";
		return "X^" + std::to_string(this->Degree);
	}
};

class MonomialBasis1D : public FunctionalBasis1D
{
private:
	int _maxPolynomialDegree;

public:
	MonomialBasis1D(int maxPolynomialDegree, CartesianGrid1D* grid, function<double(double)> sourceFunction)
		:FunctionalBasis1D(grid, sourceFunction)
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
		return "monomials_p" + std::to_string(this->_maxPolynomialDegree);
	}
};


class MonomialGlobalBasis1D : public FunctionalGlobalBasis1D
{
private:
	int _maxPolynomialDegree;

public:
	MonomialGlobalBasis1D(int maxPolynomialDegree, CartesianGrid1D* grid, function<double(double)> sourceFunction)
		:FunctionalGlobalBasis1D(grid, sourceFunction)
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
		return "globalmonomials_p" + std::to_string(this->_maxPolynomialDegree);
	}
};