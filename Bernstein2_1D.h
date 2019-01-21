#pragma once
#include "IBasisFunction.h"
#include "Utils.h"
#include <string>
#include <math.h>
#include <assert.h>
using namespace std;

class Bernstein2_1D : public IBasisFunction1D
{
private:
	int _degree;
	int _i;
	double _binomial;
public:
	static string Code() { return "bernstein2"; };

	Bernstein2_1D(int degree, int i)
	{
		this->LocalNumber = i;
		this->_degree = degree;
		this->_i = i;
		this->_binomial = Utils::Binomial(degree, i);
	}

	DefInterval DefinitionInterval() { return DefInterval::MinusOne_One(); }

	int GetDegree()
	{
		return this->_degree;
	}

	double Eval(double x)
	{
		assert(x >= -1 && x <= 1);
		return Bernstein2(x);
	}

	double EvalDerivative(double x)
	{
		assert(x >= -1 && x <= 1);
		return DBernstein2(x);
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
		//return binomial + var + "^" + to_string(this->_i) + " * (1-" + var + ")^" + to_string(this->_degree - this->_i);
	}

private:
	double Bernstein2(double x)
	{
		int n = this->_degree;
		int i = this->_i;
		return this->_binomial * pow(0.5*x + 0.5, i) * pow(0.5 - 0.5*x, n - i);
	}

	double DBernstein2(double x)
	{
		int n = this->_degree;
		int i = this->_i;
		if (n == 0)
			return 0;
		if (i == 0)
			return -this->_binomial * 0.5 * n * pow(0.5 - 0.5*x, n - 1);
		else if (n - i == 0)
			return this->_binomial * 0.5 * i * pow(0.5*x + 0.5, i - 1);
		else
			return this->_binomial * (i*0.5*pow(0.5*x + 0.5, i - 1)*pow(0.5 - 0.5*x, n - i) - (n - i)*0.5*pow(0.5*x + 0.5, i)*pow(0.5 - 0.5*x, n - i - 1));
	}
};