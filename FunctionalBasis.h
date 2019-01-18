#pragma once
#include "FunctionalBasisWithNumbers_OLD.h"
#include "FunctionalBasisWithObjects.h"
#include "BasisFunctionFactory.h"
#include "IBasisFunction.h"
#include "TensorPolynomial.h"
#include "Bernstein2D.h"
#include "Bernstein3D.h"
#include <Eigen/Sparse>

//----------//
//    1D    //
//----------//

class FunctionalBasis1DOLD : public FunctionalBasisWithNumbers
{
private:
	int _maxPolynomialDegree;
	string _basisCode;

public:
	FunctionalBasis1DOLD(string basisCode, int maxPolynomialDegree)
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

class FunctionalBasis1D : public FunctionalBasisWithObjects<IBasisFunction1D>
{
private:
	int _maxPolynomialDegree;
	string _basisCode;

public:
	FunctionalBasis1D(string basisCode, int maxPolynomialDegree)
		:FunctionalBasisWithObjects<IBasisFunction1D>()
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

	function<double(double)> GetApproximateFunction(const Eigen::VectorXd &solution, BigNumber startIndex)
	{
		function<double(double)> approximate = [this, solution, startIndex](double t) {
			double result = 0;
			for (unsigned int i = 0; i < this->_localFunctions.size(); i++)
			{
				auto phi = this->_localFunctions[i];
				result += solution(startIndex + i) * phi->Eval(t);
			}
			return result;
		};
		return approximate;
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
	bool _fullTensorization;

public:
	FunctionalBasis2D(string basisCode, int maxPolynomialDegree, bool fullTensorization)
		:FunctionalBasisWithObjects<IBasisFunction2D>()
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;
		this->_basisCode = basisCode;
		this->_fullTensorization = fullTensorization;

		int functionNumber = 0;

		if (fullTensorization)
		{
			for (int j = 0; j <= maxPolynomialDegree; j++)
			{
				IBasisFunction1D* polyY = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, j);
				for (int i = 0; i <= maxPolynomialDegree; i++)
				{
					IBasisFunction1D* polyX = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, i);
					this->_localFunctions[functionNumber] = new TensorPolynomial2D(functionNumber, polyX, polyY);
					functionNumber++;
				}
			}
		}
		else
		{
			if (basisCode.compare(Bernstein2D::Code()) == 0)
			{
				for (int j = 0; j <= maxPolynomialDegree; j++)
				{
					for (int i = 0; i <= maxPolynomialDegree - j; i++)
					{
						this->_localFunctions[functionNumber] = new Bernstein2D(functionNumber, maxPolynomialDegree, i, j);
						functionNumber++;
					}
				}
			}
			else
			{
				for (int degree = 0; degree <= maxPolynomialDegree; degree++)
				{
					for (int j = 0; j <= degree; j++)
					{
						int i = degree - j;

						IBasisFunction1D* polyX = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, i);
						IBasisFunction1D* polyY = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, j);
						this->_localFunctions[functionNumber] = new TensorPolynomial2D(functionNumber, polyX, polyY);
						functionNumber++;
					}
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
		string name = this->_basisCode + "_p" + std::to_string(this->_maxPolynomialDegree);
		if (this->_fullTensorization)
			name += "_ft";
		return name;
	}

	function<double(double, double)> GetApproximateFunction(const Eigen::VectorXd &solution, BigNumber startIndex)
	{
		function<double(double, double)> approximate = [this, solution, startIndex](double t, double u) {
			double result = 0;
			for (unsigned int i = 0; i < this->_localFunctions.size(); i++)
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
	bool _fullTensorization;

public:
	FunctionalBasis3D(string basisCode, int maxPolynomialDegree, bool fullTensorization)
		:FunctionalBasisWithObjects<IBasisFunction3D>()
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;
		this->_basisCode = basisCode;
		this->_fullTensorization = fullTensorization;

		int functionNumber = 0;

		if (fullTensorization)
		{
			for (int k = 0; k <= maxPolynomialDegree; k++)
			{
				IBasisFunction1D* polyZ = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, k);
				for (int j = 0; j <= maxPolynomialDegree; j++)
				{
					IBasisFunction1D* polyY = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, j);
					for (int i = 0; i <= maxPolynomialDegree; i++)
					{
						IBasisFunction1D* polyX = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, i);
						this->_localFunctions[functionNumber] = new TensorPolynomial3D(functionNumber, polyX, polyY, polyZ);
						functionNumber++;
					}
				}
			}
		}
		else
		{
			if (basisCode.compare(Bernstein3D::Code()) == 0)
			{
				for (int degZ = 0; degZ <= maxPolynomialDegree; degZ++)
				{
					for (int degY = 0; degY <= maxPolynomialDegree - degZ; degY++)
					{
						for (int degX = 0; degX <= maxPolynomialDegree - degY - degZ; degX++)
						{
							this->_localFunctions[functionNumber] = new Bernstein3D(functionNumber, maxPolynomialDegree, degX, degY, degZ);
							functionNumber++;
						}
					}
				}
			}
			else
			{
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
							this->_localFunctions[functionNumber] = new TensorPolynomial3D(functionNumber, polyX, polyY, polyZ);
							functionNumber++;
						}
					}
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
		string name = this->_basisCode + "_p" + std::to_string(this->_maxPolynomialDegree);
		if (this->_fullTensorization)
			name += "_ft";
		return name;
	}

	function<double(double, double, double)> GetApproximateFunction(const Eigen::VectorXd &solution, BigNumber startIndex)
	{
		function<double(double, double, double)> approximate = [this, solution, startIndex](double t, double u, double v) {
			double result = 0;
			for (unsigned int i = 0; i < this->_localFunctions.size(); i++)
			{
				auto phi = this->_localFunctions[i];
				result += solution(startIndex + i) * phi->Eval(t, u, v);
			}
			return result;
		};
		return approximate;
	}
};