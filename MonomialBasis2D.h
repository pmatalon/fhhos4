#pragma once
#include "FunctionalBasis2D.h"
#include "FunctionalGlobalBasis2D.h"
#include "IBasisFunction2D.h"
#include "IPolynomialFunction.h"
#include <math.h>

using namespace std;

class Monomial2D : public IBasisFunction2D, public IPolynomialFunction
{
public:
	int DegreeX;
	int DegreeY;

	Monomial2D(int degreeX, int degreeY)
	{
		this->DegreeX = degreeX;
		this->DegreeY = degreeY;
	}

	int GetDegree()
	{
		return this->DegreeX + this->DegreeY;
	}

	double Eval(double x, double y)
	{
		return pow(x, this->DegreeX)*pow(y, this->DegreeY);
	}

	double EvalGradX(double x, double y)
	{
		if (this->DegreeX == 0)
			return 0;
		return this->DegreeX*pow(x, this->DegreeX - 1)*pow(y, this->DegreeY);
	}

	double EvalGradY(double x, double y)
	{
		if (this->DegreeY == 0)
			return 0;
		return this->DegreeY*pow(y, this->DegreeY - 1)*pow(x, this->DegreeX);
	}

	string ToString()
	{
		if (this->DegreeX == 0 && this->DegreeY == 0)
			return "1";
		if (this->DegreeX == 1 && this->DegreeY == 0)
			return "X";
		if (this->DegreeX == 0 && this->DegreeY == 1)
			return "Y";
		if (this->DegreeX == 1 && this->DegreeY == 1)
			return "XY";
		if (this->DegreeX == 0)
			return "Y^" + std::to_string(this->DegreeY);
		if (this->DegreeY == 0)
			return "X^" + std::to_string(this->DegreeX);
		if (this->DegreeX == 1)
			return "XY^" + std::to_string(this->DegreeY);
		if (this->DegreeY == 1)
			return "X^" + std::to_string(this->DegreeX) + "Y";
		return "X^" + std::to_string(this->DegreeX) + "Y^" + std::to_string(this->DegreeY);
	}
};

class MonomialBasis2D : public FunctionalBasis2D
{
private:
	int _maxPolynomialDegree;

public:
	MonomialBasis2D(int maxPolynomialDegree, int penalizationCoefficient, function<double(double, double)> sourceFunction)
		:FunctionalBasis2D(penalizationCoefficient, sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		int functionNumber = 0;
		for (int degree = 0; degree <= this->_maxPolynomialDegree; degree++)
		{
			for (int j = 0; j <= degree; j++)
			{
				int i = degree - j;
				this->_localFunctions[functionNumber++] = new Monomial2D(i, j);
			}
		}
	}

	int GetDegree()
	{
		return this->_maxPolynomialDegree;
	}

	std::string Name()
	{
		return "monomials_p" + std::to_string(this->_maxPolynomialDegree);
	}
};


class MonomialGlobalBasis2D : public FunctionalGlobalBasis2D
{
private:
	int _maxPolynomialDegree;

public:
	MonomialGlobalBasis2D(int maxPolynomialDegree, CartesianGrid2D* grid, int penalizationCoefficient, function<double(double, double)> sourceFunction)
		:FunctionalGlobalBasis2D(grid, penalizationCoefficient, sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		int functionNumber = 0;
		for (int degree = 0; degree <= this->_maxPolynomialDegree; degree++)
		{
			for (int j = 0; j <= degree; j++)
			{
				int i = degree - j;
				this->_localFunctions[functionNumber++] = new Monomial2D(i, j);
			}
		}
	}

	int GetDegree()
	{
		return this->_maxPolynomialDegree;
	}

	std::string Name()
	{
		return "globalmonomials_p" + std::to_string(this->_maxPolynomialDegree);
	}
};