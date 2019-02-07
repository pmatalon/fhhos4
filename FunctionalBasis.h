#pragma once
#include "BasisFunctionFactory.h"
#include "BasisFunction.h"
#include "TensorPolynomial.h"
#include "Bernstein2D.h"
#include "Bernstein3D.h"
#include <Eigen/Sparse>

class FunctionalBasis
{
public:
	vector<BasisFunction*> LocalFunctions;
	virtual std::string Name() = 0;

	virtual int GetDegree() = 0;

	int NumberOfLocalFunctionsInElement(Element* element)
	{
		return static_cast<int>(this->LocalFunctions.size());
	}

	BigNumber GlobalFunctionNumber(Element* element, BasisFunction* phi)
	{
		return element->Number * static_cast<int>(this->LocalFunctions.size()) + phi->LocalNumber; // the numbers start at 0
	}

	virtual ~FunctionalBasis() {}
};

//----------//
//    1D    //
//----------//

class FunctionalBasis1D : public FunctionalBasis
{
private:
	int _maxPolynomialDegree;
	string _basisCode;

public:
	FunctionalBasis1D(string basisCode, int maxPolynomialDegree)
		:FunctionalBasis()
	{
		this->_maxPolynomialDegree = maxPolynomialDegree;
		this->_basisCode = basisCode;

		for (int i = 0; i <= maxPolynomialDegree; i++)
			this->LocalFunctions.push_back(BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, i));
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
			for (unsigned int i = 0; i < this->LocalFunctions.size(); i++)
			{
				IBasisFunction1D* phi = static_cast<IBasisFunction1D*>(this->LocalFunctions[i]);
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

class FunctionalBasis2D : public FunctionalBasis
{
private:
	int _maxPolynomialDegree;
	string _basisCode;
	bool _fullTensorization;

public:
	FunctionalBasis2D(string basisCode, int maxPolynomialDegree, bool fullTensorization)
		:FunctionalBasis()
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
					this->LocalFunctions.push_back(new TensorPolynomial2D(functionNumber, polyX, polyY));
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
						this->LocalFunctions.push_back(new Bernstein2D(functionNumber, maxPolynomialDegree, i, j));
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
						this->LocalFunctions.push_back(new TensorPolynomial2D(functionNumber, polyX, polyY));
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
			for (unsigned int i = 0; i < this->LocalFunctions.size(); i++)
			{
				IBasisFunction2D* phi = static_cast<IBasisFunction2D*>(this->LocalFunctions[i]);
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

class FunctionalBasis3D : public FunctionalBasis
{
private:
	int _maxPolynomialDegree;
	string _basisCode;
	bool _fullTensorization;

public:
	FunctionalBasis3D(string basisCode, int maxPolynomialDegree, bool fullTensorization)
		:FunctionalBasis()
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
						this->LocalFunctions.push_back(new TensorPolynomial3D(functionNumber, polyX, polyY, polyZ));
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
							this->LocalFunctions.push_back(new Bernstein3D(functionNumber, maxPolynomialDegree, degX, degY, degZ));
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
							this->LocalFunctions.push_back(new TensorPolynomial3D(functionNumber, polyX, polyY, polyZ));
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
			for (unsigned int i = 0; i < this->LocalFunctions.size(); i++)
			{
				IBasisFunction3D* phi = static_cast<IBasisFunction3D*>(this->LocalFunctions[i]);
				result += solution(startIndex + i) * phi->Eval(t, u, v);
			}
			return result;
		};
		return approximate;
	}
};