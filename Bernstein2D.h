#pragma once
#include "IBasisFunction.h"
#include "DefInterval.h"
#include "Utils.h"
#include <math.h>

class Bernstein2D : public IBasisFunction2D
{
private:
	int _degX;
	int _degY;
	int _degMixed;

	int _degree;
	int _binomial;
public:
	static string Code() { return Bernstein1D::Code(); };

	Bernstein2D(int localNumber, int totalDegree, int degX, int degY)
	{
		this->LocalNumber = localNumber;
		
		this->_degree = totalDegree;
		this->_degX = degX;
		this->_degY = degY;
		this->_degMixed = totalDegree - degX - degY;

		this->_binomial = Utils::Factorial(this->_degree) / (Utils::Factorial(this->_degX) * Utils::Factorial(this->_degY) * Utils::Factorial(this->_degMixed));
	}

	DefInterval DefinitionInterval() { return  DefInterval::MinusOne_One(); }

	int GetDegree()
	{
		return this->_degree;
	}

	double Eval(double x, double y)
	{
		// Bernstein on [-1,1]: change of variable
		return BivariateBernstein(0.5*x + 0.5, 0.5*y + 0.5);
	}

	double EvalGradX(double x, double y)
	{
		return 0.5 * GradBivariateBernsteinX(0.5*x + 0.5, 0.5*y + 0.5);
	}

	double EvalGradY(double x, double y)
	{
		return 0.5 * GradBivariateBernsteinY(0.5*x + 0.5, 0.5*y + 0.5);
	}

	string ToString()
	{
		if (this->_degree == 0)
			return "1";

		string binomialPart = this->_binomial == 1 ? "" : to_string(this->_binomial) + " * ";

		string xPart = "";
		if (this->_degX == 1)
			xPart = "X";
		else if (this->_degX > 1)
			xPart = "X^" + to_string(this->_degX);

		string yPart = "";
		if (this->_degY == 1)
			yPart = "Y";
		else if (this->_degY > 1)
			yPart = "Y^" + to_string(this->_degY);

		string mixedPart = "";
		if (this->_degMixed == 1)
			mixedPart = "(1-X-Y)";
		else if (this->_degMixed > 1)
			mixedPart = "(1-X-Y)^" + to_string(this->_degMixed);

		if (xPart.compare("") == 0)
		{
			if (yPart.compare("") == 0)
				return binomialPart + mixedPart;
			return binomialPart + (mixedPart.compare("") == 0 ? yPart : yPart + mixedPart);
		}
		if (yPart.compare("") == 0)
			return binomialPart + (mixedPart.compare("") == 0 ? xPart : xPart + mixedPart);
		if (mixedPart.compare("") == 0)
			return binomialPart + xPart + yPart;

		return binomialPart + xPart + " * " + yPart + " * " + mixedPart;
	}

private:
	// Bernstein polynomial on [0,1]
	double BivariateBernstein(double x, double y)
	{
		assert(x >= 0 && x <= 1 && y >= 0 && y <= 1);
		return this->_binomial * pow(x, this->_degX) * pow(y, this->_degY) * pow(1 - x - y, this->_degMixed);
	}

	double GradBivariateBernsteinX(double x, double y)
	{
		if (this->_degX == 0 && this->_degMixed == 0)
			return 0;

		if (this->_degX == 0)
			return -this->_binomial * pow(y, this->_degY) * this->_degMixed * pow(1 - x - y, this->_degMixed - 1);

		if (this->_degMixed == 0)
			return this->_binomial * this->_degX * pow(x, this->_degX - 1) * pow(y, this->_degY);

		return this->_binomial * pow(y, this->_degY) * (this->_degX * pow(x, this->_degX - 1) * pow(1 - x - y, this->_degMixed) - this->_degMixed * pow(1 - x - y, this->_degMixed - 1) * pow(x, this->_degX));
	}

	double GradBivariateBernsteinY(double x, double y)
	{
		if (this->_degY == 0 && this->_degMixed == 0)
			return 0;

		if (this->_degY == 0)
			return -this->_binomial * pow(x, this->_degX) * this->_degMixed * pow(1 - x - y, this->_degMixed - 1);

		if (this->_degMixed == 0)
			return this->_binomial * this->_degY * pow(x, this->_degX) * pow(y, this->_degY - 1);

		return this->_binomial * pow(x, this->_degX) * (this->_degY * pow(y, this->_degY - 1) * pow(1 - x - y, this->_degMixed) - this->_degMixed * pow(1 - x - y, this->_degMixed - 1) * pow(y, this->_degY));
	}
};