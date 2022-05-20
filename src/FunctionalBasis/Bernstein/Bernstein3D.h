#pragma once
#include "../BasisFunction.h"
#include "../../Utils/Utils.h"
#include <math.h>

class Bernstein3D : public IBasisFunction3D
{
private:
	int _degX;
	int _degY;
	int _degZ;
	int _degMixed;

	int _degree;
	int _binomial;
public:
	Bernstein3D(int localNumber, int totalDegree, int degX, int degY, int degZ)
	{
		this->LocalNumber = localNumber;

		this->_degree = totalDegree;
		this->_degX = degX;
		this->_degY = degY;
		this->_degZ = degZ;
		this->_degMixed = totalDegree - degX - degY - degZ;

		this->_binomial = Utils::Factorial(this->_degree) / (Utils::Factorial(this->_degX) * Utils::Factorial(this->_degY) * Utils::Factorial(this->_degZ) * Utils::Factorial(this->_degMixed));
	}

	int GetDegree() const
	{
		return this->_degree;
	}

	double Eval(double x, double y, double z) const
	{
		// Bernstein on [-1,1]: change of variable
		return TrivariateBernstein(0.5*x + 0.5, 0.5*y + 0.5, 0.5*z + 0.5);
	}

	double EvalGradX(double x, double y, double z) const
	{
		return 0.5 * GradTrivariateBernsteinX(0.5*x + 0.5, 0.5*y + 0.5, 0.5*z + 0.5);
	}

	double EvalGradY(double x, double y, double z) const
	{
		return 0.5 * GradTrivariateBernsteinY(0.5*x + 0.5, 0.5*y + 0.5, 0.5*z + 0.5);
	}

	double EvalGradZ(double x, double y, double z) const
	{
		return 0.5 * GradTrivariateBernsteinZ(0.5*x + 0.5, 0.5*y + 0.5, 0.5*z + 0.5);
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

		string zPart = "";
		if (this->_degZ == 1)
			zPart = "Z";
		else if (this->_degZ > 1)
			zPart = "Z^" + to_string(this->_degZ);

		string mixedPart = "";
		if (this->_degMixed == 1)
			mixedPart = "(1-X-Y-Z)";
		else if (this->_degMixed > 1)
			mixedPart = "(1-X-Y-Z)^" + to_string(this->_degMixed);

		/*if (xPart.compare("") == 0)
		{
			if (yPart.compare("") == 0)
				return binomialPart + mixedPart;
			return binomialPart + (mixedPart.compare("") == 0 ? yPart : yPart + mixedPart);
		}
		if (yPart.compare("") == 0)
			return binomialPart + (mixedPart.compare("") == 0 ? xPart : xPart + mixedPart);
		if (mixedPart.compare("") == 0)
			return binomialPart + xPart + yPart + zPart;*/

		return binomialPart + xPart + " * " + yPart + " * " + zPart + " * " + mixedPart;
	}

private:
	// Bernstein polynomial on [0,1]
	double TrivariateBernstein(double x, double y, double z) const
	{
		return this->_binomial * pow(x, this->_degX) * pow(y, this->_degY) * pow(z, this->_degZ) * pow(1 - x - y - z, this->_degMixed);
	}

	double GradTrivariateBernsteinX(double x, double y, double z) const
	{
		if (this->_degX == 0 && this->_degMixed == 0)
			return 0;

		if (this->_degX == 0)
			return -this->_binomial * pow(y, this->_degY) * pow(z, this->_degZ) * this->_degMixed * pow(1 - x - y - z, this->_degMixed - 1);

		if (this->_degMixed == 0)
			return this->_binomial * this->_degX * pow(x, this->_degX - 1) * pow(y, this->_degY) * pow(z, this->_degZ);

		return this->_binomial * pow(y, this->_degY) * pow(z, this->_degZ) * (this->_degX * pow(x, this->_degX - 1) * pow(1 - x - y - z, this->_degMixed) - this->_degMixed * pow(1 - x - y - z, this->_degMixed - 1) * pow(x, this->_degX));
	}

	double GradTrivariateBernsteinY(double x, double y, double z) const
	{
		if (this->_degY == 0 && this->_degMixed == 0)
			return 0;

		if (this->_degY == 0)
			return -this->_binomial * pow(x, this->_degX) * pow(z, this->_degZ) * this->_degMixed * pow(1 - x - y - z, this->_degMixed - 1);

		if (this->_degMixed == 0)
			return this->_binomial * this->_degY * pow(x, this->_degX) * pow(y, this->_degY - 1) * pow(z, this->_degZ);

		return this->_binomial * pow(x, this->_degX) * pow(z, this->_degZ) * (this->_degY * pow(y, this->_degY - 1) * pow(1 - x - y - z, this->_degMixed) - this->_degMixed * pow(1 - x - y - z, this->_degMixed - 1) * pow(y, this->_degY));
	}

	double GradTrivariateBernsteinZ(double x, double y, double z) const
	{
		if (this->_degZ == 0 && this->_degMixed == 0)
			return 0;

		if (this->_degZ == 0)
			return -this->_binomial * pow(x, this->_degX) * pow(y, this->_degY) * this->_degMixed * pow(1 - x - y - z, this->_degMixed - 1);

		if (this->_degMixed == 0)
			return this->_binomial * this->_degZ * pow(x, this->_degX) * pow(y, this->_degY) * pow(z, this->_degZ - 1);

		return this->_binomial * pow(x, this->_degX) * pow(y, this->_degY) * (this->_degZ * pow(z, this->_degZ - 1) * pow(1 - x - y - z, this->_degMixed) - this->_degMixed * pow(1 - x - y - z, this->_degMixed - 1) * pow(z, this->_degZ));
	}
};