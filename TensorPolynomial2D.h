#pragma once
#include "FunctionalBasis2D.h"
#include "IBasisFunction2D.h"
#include <math.h>
#include "MonomialBasis1D.h"
#include "BernsteinBasis1D.h"
#include "LegendreBasis1D.h"
using namespace std;

class TensorPolynomial2D : public IBasisFunction2D
{
private:
	IBasisFunction1D* _funcX;
	IBasisFunction1D* _funcY;
public:

	TensorPolynomial2D(IBasisFunction1D* funcX, IBasisFunction1D* funcY)
	{
		this->_funcX = funcX;
		this->_funcY = funcY;
	}

	RefInterval ReferenceInterval() { return this->_funcX->ReferenceInterval(); }

	int GetDegree()
	{
		return this->_funcX->GetDegree() + this->_funcY->GetDegree();
	}

	double Eval(double x, double y)
	{
		return this->_funcX->Eval(x) * this->_funcY->Eval(y);
	}

	double EvalGradX(double x, double y)
	{
		return this->_funcX->EvalDerivative(x) * this->_funcY->Eval(y);
	}

	double EvalGradY(double x, double y)
	{
		return this->_funcX->Eval(x) * this->_funcY->EvalDerivative(y);
	}

	string ToString()
	{
		string polyX = this->_funcX->ToString("X");
		string polyY = this->_funcY->ToString("Y");
		if (polyX.compare("1") == 0)
			return polyY;
		if (polyY.compare("1") == 0)
			return polyX;
		return polyX + " * " + polyY;
	}
};



class MonomialBasis2D : public FunctionalBasis2D
{
private:
	int _maxPolynomialDegree;

public:
	MonomialBasis2D(int maxPolynomialDegree, function<double(double, double)> sourceFunction)
		:FunctionalBasis2D(sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		int functionNumber = 0;
		for (int degree = 0; degree <= this->_maxPolynomialDegree; degree++)
		{
			for (int j = 0; j <= degree; j++)
			{
				int i = degree - j;

				IBasisFunction1D* polyX = new Monomial1D(i);
				IBasisFunction1D* polyY = new Monomial1D(j);
				this->_localFunctions[functionNumber++] = new TensorPolynomial2D(polyX, polyY);
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



class LegendreBasis2D : public FunctionalBasis2D
{
private:
	int _maxPolynomialDegree;

public:
	LegendreBasis2D(int maxPolynomialDegree, function<double(double, double)> sourceFunction)
		:FunctionalBasis2D(sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		int functionNumber = 0;
		for (int degree = 0; degree <= this->_maxPolynomialDegree; degree++)
		{
			for (int j = 0; j <= degree; j++)
			{
				int i = degree - j;

				IBasisFunction1D* polyX = new Legendre1D(i, false);
				IBasisFunction1D* polyY = new Legendre1D(j, false);
				this->_localFunctions[functionNumber++] = new TensorPolynomial2D(polyX, polyY);
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


class BernsteinBasis2D : public FunctionalBasis2D
{
private:
	int _maxPolynomialDegree;

public:
	BernsteinBasis2D(int maxPolynomialDegree, function<double(double, double)> sourceFunction)
		:FunctionalBasis2D(sourceFunction)
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;

		int functionNumber = 0;
		for (int degree = 0; degree <= this->_maxPolynomialDegree; degree++)
		{
			for (int j = 0; j <= degree; j++)
			{
				int i = degree - j;

				IBasisFunction1D* polyX = new Bernstein1D(_maxPolynomialDegree, i);
				IBasisFunction1D* polyY = new Bernstein1D(_maxPolynomialDegree, j);
				this->_localFunctions[functionNumber++] = new TensorPolynomial2D(polyX, polyY);
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
};