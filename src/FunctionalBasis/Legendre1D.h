#pragma once
#include "BasisFunction.h"
#include <string>
#include <math.h>
#include <assert.h>
using namespace std;

class Legendre1D : public IBasisFunction1D
{
public:
	int Degree;

	static string Code() { return "legendre"; };

	Legendre1D(int degree)
	{
		this->LocalNumber = degree;
		this->Degree = degree;
	}

	int GetDegree() const
	{
		return this->Degree;
	}

	virtual double Eval(double x) override
	{
		this->TestIsInReferenceInterval(x);
		return Legendre(this->Degree, x);
	}

	virtual double EvalDerivative(double x) override
	{
		this->TestIsInReferenceInterval(x);
		return DLegendre(this->Degree, x);
	}

	static double Legendre(int n, double x)
	{
		if (n == 0)
			return 1;
		if (n == 1)
			return x;
		// Recurrence (Bonnet's formula)
		return ((2 * n - 1)*x*Legendre(n - 1, x) - (n - 1)*Legendre(n - 2, x))/n;
	}

	static double DLegendre(int n, double x)
	{
		if (n == 0)
			return 0;
		if (n == 1)
			return 1;
		/*if (n == 2)
			return 3 * x;
		if (n == 3)
			return 7.5 * pow(x, 2) - 1.5;
		if (n == 4)
			return 17.5 * pow(x, 3) - 7.5 * x;*/
		//dy(i + 2) = ((2 * i + 1)*y(i + 1) + (2 * i + 1)*x*dy(i + 1) - i * dy(i)) / (i + 1);

		//return inverseNorm*DLegendre(n - 2, x) + (2 * n - 1)*Legendre(n - 1, x);
		return ((2 * n - 1)*Legendre(n - 1, x) + (2 * n - 1)*x*DLegendre(n - 1, x) - (n-1) * DLegendre(n - 2, x)) / n;
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
		return "Legendre(" + std::to_string(this->Degree) + ", " + var + ")";
	}
};

class NormalizedLegendre1D : public Legendre1D
{
private:
	double _inverseNorm;
public:
	static string Code() { return "nlegendre"; };

	NormalizedLegendre1D(int degree) : Legendre1D(degree)
	{
		_inverseNorm = sqrt(this->Degree + 0.5);
	}

	virtual double Eval(double x) override
	{
		return _inverseNorm * Legendre1D::Eval(x);
	}

	virtual double EvalDerivative(double x) override
	{
		return _inverseNorm * Legendre1D::EvalDerivative(x);
	}
};