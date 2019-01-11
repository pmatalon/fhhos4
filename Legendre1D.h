#pragma once
#include "IBasisFunction.h"
#include <string>
#include <math.h>
#include <assert.h>
using namespace std;

class Legendre1D : public IBasisFunction1D
{
public:
	int Degree;
	bool Normalized;

	static string Code() { return "legendre"; };

	Legendre1D(int degree)
	{
		this->Degree = degree;
		this->Normalized = false;
	}

	Legendre1D(int degree, bool normalized)
	{
		this->Degree = degree;
		this->Normalized = normalized;
	}

	RefInterval ReferenceInterval() { return RefInterval::MinusOne_One(); }

	int GetDegree()
	{
		return this->Degree;
	}

	double Eval(double x)
	{
		assert(x >= -1 && x <= 1);
		double value = Legendre(this->Degree, x);
		if (this->Normalized)
			return sqrt(this->Degree + 0.5) * value;
		return value;
	}

	double EvalDerivative(double x)
	{
		assert(x >= -1 && x <= 1);
		double value = DLegendre(this->Degree, x);
		if (this->Normalized)
			return sqrt(this->Degree + 0.5) * value;
		return value;
	}

	static double Legendre(int n, double x)
	{
		if (n == 0)
			return 1;
		if (n == 1)
			return x;
		/*if (n == 2)
			return 1.5 * pow(x, 2) - 0.5;
		if (n == 3)
			return 2.5 * pow(x, 3) - 1.5*x;
		if (n == 4)
			return 0.125 * (35 * pow(x, 4) - 30 * pow(x, 2) + 3);*/
		
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