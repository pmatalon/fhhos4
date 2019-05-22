#pragma once
#include "BasisFunctionFactory.h"
#include "BasisFunction.h"
#include "TensorPolynomial.h"
#include "Bernstein2D.h"
#include "Bernstein3D.h"
#include <Eigen/Sparse>
template <int Dim>
class Element;

template <int Dim>
class FunctionalBasis
{
protected:
	int _maxPolynomialDegree;
	string _basisCode;

public:
	bool FullTensorization;
	vector<BasisFunction<Dim>*> LocalFunctions;

	FunctionalBasis(string basisCode, int maxPolynomialDegree) 
		: FunctionalBasis(basisCode, maxPolynomialDegree, false)
	{}

	FunctionalBasis(string basisCode, int maxPolynomialDegree, bool fullTensorization)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);
		this->_maxPolynomialDegree = maxPolynomialDegree;
		this->_basisCode = basisCode;
		this->FullTensorization = fullTensorization;

		//----------//
		//    1D    //
		//----------//
		if (Dim == 1)
		{
			this->FullTensorization = false;
			for (int i = 0; i <= maxPolynomialDegree; i++)
			{
				BasisFunction<Dim>* poly = dynamic_cast<BasisFunction<Dim>*>(BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, i));
				this->LocalFunctions.push_back(poly);
			}
		}
		//----------//
		//    2D    //
		//----------//
		else if (Dim == 2)
		{
			int functionNumber = 0;

			if (fullTensorization)
			{
				for (int j = 0; j <= maxPolynomialDegree; j++)
				{
					for (int i = 0; i <= maxPolynomialDegree; i++)
					{
						IBasisFunction1D* polyX = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, i);
						IBasisFunction1D* polyY = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, j);
						BasisFunction<Dim>* poly2D = dynamic_cast<BasisFunction<Dim>*>(new TensorPolynomial2D(functionNumber, polyX, polyY));
						this->LocalFunctions.push_back(poly2D);
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
							BasisFunction<Dim>* poly2D = dynamic_cast<BasisFunction<Dim>*>(new Bernstein2D(functionNumber, maxPolynomialDegree, i, j));
							this->LocalFunctions.push_back(poly2D);
							functionNumber++;
						}
					}
				}
				else if (basisCode.compare(Hemker1D::Code()) == 0)
				{
					for (int degree = 0; degree <= maxPolynomialDegree; degree++)
					{
						for (int j = 0; j <= degree; j++)
						{
							int i = degree - j;

							IBasisFunction1D* polyX = BasisFunctionFactory::Create(basisCode, i, i);
							IBasisFunction1D* polyY = BasisFunctionFactory::Create(basisCode, j, j);

							BasisFunction<Dim>* poly2D = dynamic_cast<BasisFunction<Dim>*>(new TensorPolynomial2D(functionNumber, polyX, polyY));
							this->LocalFunctions.push_back(poly2D);
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

							BasisFunction<Dim>* poly2D = dynamic_cast<BasisFunction<Dim>*>(new TensorPolynomial2D(functionNumber, polyX, polyY));
							this->LocalFunctions.push_back(poly2D);
							functionNumber++;
						}
					}
				}
			}
		}
		//----------//
		//    3D    //
		//----------//
		else if (Dim == 3)
		{
			int functionNumber = 0;

			if (fullTensorization)
			{
				for (int k = 0; k <= maxPolynomialDegree; k++)
				{
					for (int j = 0; j <= maxPolynomialDegree; j++)
					{
						for (int i = 0; i <= maxPolynomialDegree; i++)
						{
							IBasisFunction1D* polyX = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, i);
							IBasisFunction1D* polyY = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, j);
							IBasisFunction1D* polyZ = BasisFunctionFactory::Create(basisCode, maxPolynomialDegree, k);

							BasisFunction<Dim>* poly3D = dynamic_cast<BasisFunction<Dim>*>(new TensorPolynomial3D(functionNumber, polyX, polyY, polyZ));
							this->LocalFunctions.push_back(poly3D);
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
								BasisFunction<Dim>* poly3D = dynamic_cast<BasisFunction<Dim>*>(new Bernstein3D(functionNumber, maxPolynomialDegree, degX, degY, degZ));
								this->LocalFunctions.push_back(poly3D);
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

								BasisFunction<Dim>* poly3D = dynamic_cast<BasisFunction<Dim>*>(new TensorPolynomial3D(functionNumber, polyX, polyY, polyZ));
								this->LocalFunctions.push_back(poly3D);
								functionNumber++;
							}
						}
					}
				}
			}
		}
	}

	string Name()
	{
		string name = this->_basisCode + "_p" + std::to_string(this->_maxPolynomialDegree);
		if (this->FullTensorization)
			name += "_ft";
		return name;
	}

	int GetDegree()
	{
		return this->_maxPolynomialDegree;
	}

	int NumberOfLocalFunctionsInElement(Element<Dim>* element)
	{
		return static_cast<int>(this->LocalFunctions.size());
	}

	int Size()
	{
		return NumberOfLocalFunctionsInElement(NULL);
	}

	BigNumber GlobalFunctionNumber(Element<Dim>* element, BasisFunction<Dim>* phi)
	{
		return element->Number * static_cast<int>(this->LocalFunctions.size()) + phi->LocalNumber; // the numbers start at 0
	}

	function<double(RefPoint)> GetApproximateFunction(const Eigen::VectorXd &solution, BigNumber startIndex)
	{
		function<double(RefPoint)> approximate = [this, solution, startIndex](RefPoint p) {
			double result = 0;
			for (unsigned int i = 0; i < this->LocalFunctions.size(); i++)
			{
				BasisFunction<Dim>* phi = this->LocalFunctions[i];
				result += solution(startIndex + i) * phi->Eval(p);
			}
			return result;
		};
		return approximate;
	}

	virtual ~FunctionalBasis() 
	{
		for (auto phi : LocalFunctions)
			delete phi;
	}
};