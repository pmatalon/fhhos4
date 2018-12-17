#pragma once
#include "FunctionalBasis1D.h"
#include "IBasisFunction1D.h"
#include "IPolynomialFunction.h"
#include "Utils.h"
#include <math.h>
using namespace std;

class Bernstein1D : public IBasisFunction1D, public IPolynomialFunction
{
private:
	int _degree;
	int _i;
	double _binomial;
public:

	Bernstein1D(int degree, int i)
	{
		this->_degree = degree;
		this->_i = i;
		this->_binomial = Utils::Binomial(degree, i);
	}

	RefInterval ReferenceInterval() { return RefInterval::Zero_One(); }

	int GetDegree()
	{
		return this->_degree;
	}

	double Eval(double x)
	{
		double value = Bernstein(x);
		return value;
	}

	double EvalDerivative(double x)
	{
		double value = DBernstein(x);
		return value;
	}

	string ToString()
	{
		return "binomial(" + std::to_string(this->_degree) + ", " + std::to_string(this->_i) + ") * X^" + std::to_string(this->_i) + " * (1-X)^" + std::to_string(this->_degree - this->_i);
	}

private:
	double Bernstein(double x)
	{
		int n = this->_degree;
		int i = this->_i;
		return this->_binomial * pow(x, i) * pow(1 - x, n - i);
	}

	double DBernstein(double x)
	{
		int n = this->_degree;
		int i = this->_i;
		if (n == 0)
			return 0;
		if (i == 0)
			return -this->_binomial * n * pow(1 - x, n - 1);
		else if (n - i == 0)
			return this->_binomial * i * pow(x, i - 1);
		else
			return this->_binomial * (i*pow(x, i - 1)*pow(1 - x, n - i) - (n-i)*pow(x, i)*pow(1 - x, n - i - 1));
	}
};

class BernsteinBasis1D : public FunctionalBasis1D
{
private:
	int _maxPolynomialDegree;

public:
	BernsteinBasis1D(int maxPolynomialDegree, CartesianGrid1D* grid, function<double(double)> sourceFunction)
		:FunctionalBasis1D(grid, sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		for (int i = 0; i <= maxPolynomialDegree; i++)
			this->_localFunctions[i] = new Bernstein1D(maxPolynomialDegree, i);
	}

	int GetDegree()
	{
		return this->_maxPolynomialDegree;
	}

	string Name()
	{
		return "bernstein_p" + std::to_string(this->_maxPolynomialDegree);
	}
};