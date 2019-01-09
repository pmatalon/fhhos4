#pragma once
#include "FunctionalBasisWithNumbers.h"
#include "FunctionalBasisWithObjects.h"
#include "BasisFunctionFactory.h"
#include "IBasisFunction.h"
#include "TensorPolynomial.h"
#include <Eigen/Sparse>

//----------//
//    1D    //
//----------//

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

//----------//
//    2D    //
//----------//

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

	function<double(double, double)> GetApproximateFunction(const Eigen::VectorXd &solution, BigNumber startIndex)
	{
		function<double(double, double)> approximate = [this, solution, startIndex](double t, double u) {
			double result = 0;
			for (int i = 0; i < this->_localFunctions.size(); i++)
			{
				auto phi = this->_localFunctions[i];
				result += solution(startIndex + i) * phi->Eval(t, u);
			}
			return result;
		};
		return approximate;
	}
};

//----------//
//    3D    //
//----------//

class FunctionalBasis3D : public FunctionalBasisWithObjects<IBasisFunction3D>
{
private:
	int _maxPolynomialDegree;
	string _basisCode;

public:
	FunctionalBasis3D(string basisCode, int maxPolynomialDegree)
		:FunctionalBasisWithObjects<IBasisFunction3D>()
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;
		this->_basisCode = basisCode;

		int functionNumber = 0;
		for (int degree = 0; degree <= maxPolynomialDegree; degree++)
		{
			for (int degZ = 0; degZ <= degree; degZ++)
			{
				for (int degY = 0; degY <= degree - degZ; degY++)
				{
					int degX = degree - degZ - degY;

					IBasisFunction1D* polyX = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, degX);
					IBasisFunction1D* polyY = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, degY);
					IBasisFunction1D* polyZ = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, degZ);
					this->_localFunctions[functionNumber++] = new TensorPolynomial3D(polyX, polyY, polyZ);
				}
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

	function<double(double, double, double)> GetApproximateFunction(const Eigen::VectorXd &solution, BigNumber startIndex)
	{
		function<double(double, double, double)> approximate = [this, solution, startIndex](double t, double u, double v) {
			double result = 0;
			for (int i = 0; i < this->_localFunctions.size(); i++)
			{
				auto phi = this->_localFunctions[i];
				result += solution(startIndex + i) * phi->Eval(t, u, v);
			}
			return result;
		};
		return approximate;
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