#pragma once
#include "FunctionalBasis2D.h"
#include "IBasisFunction2D.h"
#include <math.h>
#include "LegendreBasis1D.h"
using namespace std;

class Legendre2D : public IBasisFunction2D
{
public:
	int DegreeX;
	int DegreeY;
	bool Normalized;

	Legendre2D(int degreeX, int degreeY)
	{
		this->DegreeX = degreeX;
		this->DegreeY = degreeY;
	}

	Legendre2D(int degreeX, int degreeY, bool normalized)
	{
		this->DegreeX = degreeX;
		this->DegreeY = degreeY;
		this->Normalized = normalized;
	}

	RefInterval ReferenceInterval() { return RefInterval::MinusOne_One(); }

	int GetDegree()
	{
		return this->DegreeX + this->DegreeY;
	}

	double Eval(double x, double y)
	{
		double value = Legendre1D::Legendre(this->DegreeX, x) * Legendre1D::Legendre(this->DegreeY, y);
		if (this->Normalized)
			return sqrt(this->DegreeX + 0.5) * value;
		return value;
	}

	double EvalGradX(double x, double y)
	{
		double value = Legendre1D::DLegendre(this->DegreeX, x) * Legendre1D::Legendre(this->DegreeY, y);
		if (this->Normalized)
			return sqrt(this->DegreeX + 0.5) * value;
		return value;
	}

	double EvalGradY(double x, double y)
	{
		double value = Legendre1D::Legendre(this->DegreeX, x) * Legendre1D::DLegendre(this->DegreeY, y);
		if (this->Normalized)
			return sqrt(this->DegreeX + 0.5) * value;
		return value;
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
		return "Legendre(" + std::to_string(this->DegreeX) + ", X) * Legendre(" + std::to_string(this->DegreeY) + ", Y)";
	}
};

class LegendreBasis2D : public FunctionalBasis2D
{
private:
	int _maxPolynomialDegree;

public:
	LegendreBasis2D(int maxPolynomialDegree, int penalizationCoefficient, function<double(double, double)> sourceFunction)
		:FunctionalBasis2D(penalizationCoefficient, sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		int functionNumber = 0;
		for (int degree = 0; degree <= this->_maxPolynomialDegree; degree++)
		{
			for (int j = 0; j <= degree; j++)
			{
				int i = degree - j;
				this->_localFunctions[functionNumber++] = new Legendre2D(i, j);
			}
		}
	}

	int GetDegree()
	{
		return this->_maxPolynomialDegree;
	}

	string Name()
	{
		return "legendre_p" + std::to_string(this->_maxPolynomialDegree);
	}
};