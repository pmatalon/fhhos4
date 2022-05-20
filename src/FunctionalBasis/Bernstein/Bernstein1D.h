#pragma once
#include "../BasisFunction.h"
#include "../../Utils/Utils.h"
#include <math.h>
#include <assert.h>
using namespace std;

class Bernstein1D : public IBasisFunction1D
{
private:
	int _degree;
	int _i;
	double _binomial;
public:
	Bernstein1D(int degree, int i)
	{
		this->LocalNumber = i;
		this->_degree = degree;
		this->_i = i;
		this->_binomial = Utils::Binomial(degree, i);
	}

	int GetDegree() const
	{
		return this->_degree;
	}

	double Eval(double x) const
	{
		this->TestIsInReferenceInterval(x);
		// Bernstein on [-1,1]: change of variable
		return Bernstein(0.5*x + 0.5);
	}

	double EvalDerivative(double x) const
	{
		this->TestIsInReferenceInterval(x);
		return 0.5 * DBernstein(0.5*x + 0.5);
	}

	string ToString()
	{
		return this->ToString("X");
	}

	string ToString(string var)
	{
		string binomial = "";
		if (this->_binomial != 1)
			binomial = to_string(static_cast<int>(this->_binomial)) + " * ";
		string firstTerm = "";
		if (this->_i != 0)
			firstTerm = var + "^" + to_string(this->_i);
		string secondTerm = "";
		if (this->_degree - this->_i != 0)
			secondTerm = "(1-" + var + ")^" + to_string(this->_degree - this->_i);
		if (firstTerm.empty())
			return binomial + secondTerm;
		if (secondTerm.empty())
			return binomial + firstTerm;
		return binomial + firstTerm + " * " + secondTerm;
	}

private:
	// Bernstein polynomial on [0,1]
	double Bernstein(double x) const
	{
		int n = this->_degree;
		int i = this->_i;
		return this->_binomial * pow(x, i) * pow(1 - x, n - i);
	}

	double DBernstein(double x) const
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
			return this->_binomial * (i*pow(x, i - 1)*pow(1 - x, n - i) - (n - i)*pow(x, i)*pow(1 - x, n - i - 1));
	}
};