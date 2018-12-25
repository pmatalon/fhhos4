/*#pragma once
#include "FunctionalBasis2D.h"
#include "IBasisFunction2D.h"
#include <math.h>
#include "BernsteinBasis1D.h"
using namespace std;

class Bernstein2D : public IBasisFunction2D
{
private:
	Bernstein1D* _bernsteinX;
	Bernstein1D* _bernsteinY;
public:

	Bernstein2D(int degreeX, int i, int degreeY, int j)
	{
		this->_bernsteinX = new Bernstein1D(degreeX, i);
		this->_bernsteinY = new Bernstein1D(degreeY, j);
	}

	RefInterval ReferenceInterval() { return RefInterval::Zero_One(); }

	int GetDegree()
	{
		return this->_bernsteinX->GetDegree() + this->_bernsteinY->GetDegree();
	}

	double Eval(double x, double y)
	{
		return this->_bernsteinX->Eval(x) * this->_bernsteinY->Eval(y);
	}

	double EvalGradX(double x, double y)
	{
		return this->_bernsteinX->EvalDerivative(x) * this->_bernsteinY->Eval(y);
	}

	double EvalGradY(double x, double y)
	{
		return this->_bernsteinX->Eval(x) * this->_bernsteinY->EvalDerivative(y);
	}

	string ToString()
	{
		return NULL;
	}
};

class BernsteinBasis2D : public FunctionalBasis2D
{
private:
	int _maxPolynomialDegree;

public:
	BernsteinBasis2D(int maxPolynomialDegree, int penalizationCoefficient, function<double(double, double)> sourceFunction)
		:FunctionalBasis2D(penalizationCoefficient, sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		int functionNumber = 0;
		for (int degree = 0; degree <= this->_maxPolynomialDegree; degree++)
		{
			for (int j = 0; j <= degree; j++)
			{
				int i = degree - j;
				this->_localFunctions[functionNumber++] = new Bernstein2D(_maxPolynomialDegree, i, _maxPolynomialDegree, j);
			}
		}
	}

	int GetDegree()
	{
		return this->_maxPolynomialDegree;
	}

	string Name()
	{
		return "bernstein_p" + std::to_string(this->_maxPolynomialDegree);
	}
};*/