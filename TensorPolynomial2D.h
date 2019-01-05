#pragma once
#include "IBasisFunction2D.h"
#include <math.h>
#include "MonomialBasis1D.h"
#include "BernsteinBasis1D.h"
#include "Bernstein2Basis1D.h"
#include "LegendreBasis1D.h"
#include "BasisFunctionFactory.h"
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

class TensorPolynomialBasis2D : public FunctionalBasisWithObjects<IBasisFunction2D>
{
private:
	int _maxPolynomialDegree;
	string _basisCode;

public:
	TensorPolynomialBasis2D(string basisCode, int maxPolynomialDegree)
		:FunctionalBasisWithObjects<IBasisFunction2D>()
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;
		this->_basisCode = basisCode;

		int functionNumber = 0;
		for (int degree = 0; degree <= maxPolynomialDegree; degree++)
		{
			for (int j = 0; j <= degree; j++)
			{
				int i = degree - j;

				IBasisFunction1D* polyX = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, i);
				IBasisFunction1D* polyY = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, j);
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
		return this->_basisCode + "_p" + std::to_string(this->_maxPolynomialDegree);
	}
};

class MonomialBasis2D : public TensorPolynomialBasis2D
{
public:
	MonomialBasis2D(int maxPolynomialDegree)
		:TensorPolynomialBasis2D(Monomial1D::Code(), maxPolynomialDegree) {}
};
class LegendreBasis2D : public TensorPolynomialBasis2D
{
public:
	LegendreBasis2D(int maxPolynomialDegree)
		:TensorPolynomialBasis2D(Legendre1D::Code(), maxPolynomialDegree) {}
};
class BernsteinBasis2D : public TensorPolynomialBasis2D
{
public:
	BernsteinBasis2D(int maxPolynomialDegree)
		:TensorPolynomialBasis2D(Bernstein1D::Code(), maxPolynomialDegree) {}
};
class Bernstein2Basis2D : public TensorPolynomialBasis2D
{
public:
	Bernstein2Basis2D(int maxPolynomialDegree)
		:TensorPolynomialBasis2D(Bernstein2_1D::Code(), maxPolynomialDegree) {}
};