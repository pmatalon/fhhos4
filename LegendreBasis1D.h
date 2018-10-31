#pragma once
#include "FunctionalBasis1D.h"
#include "BasisFunction1D.h"
#include <math.h>
using namespace std;

class Legendre1D : public BasisFunction1D
{
public:
	int Degree;

	Legendre1D(int degree)
	{
		this->Degree = degree;
	}

	double Eval(double x)
	{
		return Legendre(this->Degree, x);
	}

	double EvalGrad(double x)
	{
		return DLegendre(this->Degree, x);
	}

private:
	static double Legendre(int n, double x)
	{
		if (n == 0)
			return 1;
		if (n == 1)
			return x;
		if (n == 2)
			return 1.5 * pow(x, 2) - 0.5;
		if (n == 3)
			return 2.5 * pow(x, 3) - 1.5*x;
		if (n == 4)
			return 0.125 * (35 * pow(x, 4) - 30 * pow(x, 2) + 3);

		return ((2 * (double)n - 1)*x*Legendre(n - 1, x) - ((double)n - 1)*Legendre(n - 2, x))/n;
	}

	static double DLegendre(int n, double x)
	{
		if (n == 0)
			return 0;
		if (n == 1)
			return 1;
		if (n == 2)
			return 3 * x;
		if (n == 3)
			return 7.5 * pow(x, 2) - 1.5;
		if (n == 4)
			return 17.5 * pow(x, 3) - 7.5 * x;

		return DLegendre(n - 2, x) + (2 * (double)n - 1)*Legendre(n - 1, x);
	}
};

class LegendreBasis1D : public FunctionalBasis1D
{
private:
	int _maxPolynomialDegree;

public:
	LegendreBasis1D(int maxPolynomialDegree, CartesianGrid1D* grid, int penalizationCoefficient, function<double(double)> sourceFunction)
		:FunctionalBasis1D(grid, penalizationCoefficient, sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		for (int i = 0; i <= maxPolynomialDegree; i++)
			this->_localFunctions[i] = new Legendre1D(i);
	}

	string Name()
	{
		return "legendre_p" + std::to_string(this->_maxPolynomialDegree);
	}
};