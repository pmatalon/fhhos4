#pragma once
#include "FunctionalBasisWithNumbers.h"
#include "FunctionalBasisWithObjects.h"
#include "BasisFunctionFactory.h"
#include "IBasisFunction.h"
#include "TensorPolynomial2D.h"
//#include "MonomialBasis1D.h"
//#include "BernsteinBasis1D.h"
//#include "Bernstein2Basis1D.h"
//#include "LegendreBasis1D.h"

class FunctionalBasis1D : public FunctionalBasisWithNumbers
{
private:
	int _maxPolynomialDegree;
	string _basisCode;

public:
	FunctionalBasis1D(string basisCode, int maxPolynomialDegree)
		:FunctionalBasisWithNumbers()
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;
		this->_basisCode = basisCode;

		for (int i = 0; i <= maxPolynomialDegree; i++)
			this->_localFunctions[i] = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, i);
	}

	int GetDegree()
	{
		return this->_maxPolynomialDegree;
	}

	string Name()
	{
		return this->_basisCode + "_p" + std::to_string(this->_maxPolynomialDegree);
	}
};

class FunctionalBasis2D : public FunctionalBasisWithObjects<IBasisFunction2D>
{
private:
	int _maxPolynomialDegree;
	string _basisCode;

public:
	FunctionalBasis2D(string basisCode, int maxPolynomialDegree)
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

/*// Monomial //
class MonomialBasis1D : public FunctionalBasis1D
{
public:
	MonomialBasis1D(int maxPolynomialDegree)
		:FunctionalBasis1D(Monomial1D::Code(), maxPolynomialDegree) {}
};
class MonomialBasis2D : public FunctionalBasis2D
{
public:
	MonomialBasis2D(int maxPolynomialDegree)
		:FunctionalBasis2D(Monomial1D::Code(), maxPolynomialDegree) {}
};

// Legendre //
class LegendreBasis1D : public FunctionalBasis1D
{
public:
	LegendreBasis1D(int maxPolynomialDegree)
		:FunctionalBasis1D(Legendre1D::Code(), maxPolynomialDegree) {}
};
class LegendreBasis2D : public FunctionalBasis2D
{
public:
	LegendreBasis2D(int maxPolynomialDegree)
		:FunctionalBasis2D(Legendre1D::Code(), maxPolynomialDegree) {}
};

// Bernstein //
class BernsteinBasis1D : public FunctionalBasis1D
{
public:
	BernsteinBasis1D(int maxPolynomialDegree)
		:FunctionalBasis1D(Bernstein1D::Code(), maxPolynomialDegree) {}
};
class BernsteinBasis2D : public FunctionalBasis2D
{
public:
	BernsteinBasis2D(int maxPolynomialDegree)
		:FunctionalBasis2D(Bernstein1D::Code(), maxPolynomialDegree) {}
};

// Bernstein 2 //
class Bernstein2Basis1D : public FunctionalBasis1D
{
public:
	Bernstein2Basis1D(int maxPolynomialDegree)
		:FunctionalBasis1D(Bernstein2_1D::Code(), maxPolynomialDegree) {}
};
class Bernstein2Basis2D : public FunctionalBasis2D
{
public:
	Bernstein2Basis2D(int maxPolynomialDegree)
		:FunctionalBasis2D(Bernstein2_1D::Code(), maxPolynomialDegree) {}
};*/