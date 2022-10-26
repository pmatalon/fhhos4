#pragma once
#include "BasisFunction.h"

template <int Dim>
class FunctionalBasis
{
protected:
	FunctionalBasis() {}
public:
	/*FunctionalBasis(int maxPolynomialDegree, bool usePolynomialSpaceQ)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);
		this->_maxPolynomialDegree = maxPolynomialDegree;
		this->UsePolynomialSpaceQ = usePolynomialSpaceQ;
		//this->IsHierarchical = basisCode.compare(Monomial1D::Code()) == 0 || basisCode.compare(Legendre1D::Code()) == 0;
		//if (basisCode.compare(Legendre1D::Code()) == 0)
			//this->IsOrthogonalOnCartesianShapes = true;

		//----------//
		//    0D    //
		//----------//
		if (Dim == 0)
		{
			this->_maxPolynomialDegree = 0;
			this->UsePolynomialSpaceQ = false;
			BasisFunction<Dim>* poly = dynamic_cast<BasisFunction<Dim>*>(new BasisFunction0D());
			this->LocalFunctions.push_back(poly);
		}
		//----------//
		//    1D    //
		//----------//
		if (Dim == 1)
		{
			this->UsePolynomialSpaceQ = false;
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

			if (usePolynomialSpaceQ)
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
				if (basisCode.compare("lagrange") == 0)
				{
					if (maxPolynomialDegree != 1)
						Utils::FatalError("Only piecewise linear approximation is implemented with the Lagrange basis.");

					BasisFunction<Dim>* lagrangeNode1 = dynamic_cast<BasisFunction<Dim>*>(new LagrangeP1_Node1());
					BasisFunction<Dim>* lagrangeNode2 = dynamic_cast<BasisFunction<Dim>*>(new LagrangeP1_Node2());
					BasisFunction<Dim>* lagrangeNode3 = dynamic_cast<BasisFunction<Dim>*>(new LagrangeP1_Node3());
					this->LocalFunctions.push_back(lagrangeNode1);
					this->LocalFunctions.push_back(lagrangeNode2);
					this->LocalFunctions.push_back(lagrangeNode3);
				}
				else if (basisCode.compare(Bernstein2D::Code()) == 0)
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
#ifdef ENABLE_3D
		else if (Dim == 3)
		{
			int functionNumber = 0;

			if (usePolynomialSpaceQ)
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
#endif // ENABLE_3D
	}*/

	virtual vector<BasisFunction<Dim>*> LocalFunctions() = 0;
	virtual string BasisCode()      const = 0;
	virtual int    GetDegree()      const = 0;
	virtual bool   IsHierarchical() const = 0;
	virtual int    Size()           const = 0;
	virtual FunctionalBasis<Dim>* CreateSameBasisForDegree(int degree) = 0;
	virtual FunctionalBasis<Dim>* CreateLowerDegreeBasis  (int degree) = 0;

	virtual bool UsePolynomialSpaceQ() const 
	{
		if (GetDegree() == 0) return false;
		return pow(GetDegree() + 1, Dim) == Size(); 
	}
	virtual bool IsOrthogonalOnCartesianShapes() const { return false; }

	string Name() const
	{
		string name = BasisCode() + "_p" + std::to_string(GetDegree());
		if (UsePolynomialSpaceQ())
			name += "_Q";
		return name;
	}

	RefFunction GetApproximateFunction(const Vector& solution, BigNumber startIndex)
	{
		RefFunction approximate = [this, &solution, startIndex](const RefPoint& p) {
			double result = 0;
			for (BasisFunction<Dim>* phi : this->LocalFunctions())
				result += solution(startIndex + phi->LocalNumber) * phi->Eval(p);
			return result;
		};
		return approximate;
	}

	virtual ~FunctionalBasis() 
	{}
};

class FunctionalBasis0D : public FunctionalBasis<0>
{
private:
	BasisFunction0D _localFunction;
public:
	vector<BasisFunction<0>*> LocalFunctions() override { return { &_localFunction }; }
	string BasisCode()           const override { return "basis0D"; }
	int    GetDegree()           const override { return 0; }
	bool   IsHierarchical()      const override { return true; }
	int    Size()                const override { return 1; }
	bool   UsePolynomialSpaceQ() const override { return false; }
	FunctionalBasis<0>* CreateSameBasisForDegree(int degree) override { assert(false); return nullptr; }
	FunctionalBasis<0>* CreateLowerDegreeBasis  (int degree) override { assert(false); return nullptr; }
};