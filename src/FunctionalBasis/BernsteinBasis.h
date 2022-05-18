#pragma once
#include "FunctionalBasis.h"
#include "TensorPolynomial.h"
#include "Bernstein1D.h"
#include "Bernstein2D.h"
#ifdef ENABLE_3D
#include "Bernstein3D.h"
#endif

template <int Dim>
class BernsteinBasis : public FunctionalBasis<Dim>
{
public:
	bool   IsHierarchical()                const override { return false; }
	string BasisCode()                     const override { return "bernstein"; }
	static string Code()                                  { return "bernstein"; };
};

class BernsteinBasis1D : public BernsteinBasis<1>
{
private:
	vector<Bernstein1D> _localFunctions;
private:
	BernsteinBasis1D() {}
public:
	BernsteinBasis1D(int maxPolynomialDegree)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);

		_localFunctions.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_localFunctions.emplace_back(maxPolynomialDegree, i);

	}
	vector<BasisFunction<1>*> LocalFunctions() override
	{
		vector<BasisFunction<1>*> list;
		list.reserve(_localFunctions.size());
		for (Bernstein1D& m : _localFunctions)
			list.push_back(&m);
		return list;
	}

	int GetDegree() const override
	{
		return Size() - 1;
	}

	int Size() const override
	{
		return (int)_localFunctions.size();
	}

	bool UsePolynomialSpaceQ() const override
	{
		return false;
	}

	FunctionalBasis<1>* CreateSameBasisForDegree(int degree) override
	{
		return new BernsteinBasis1D(degree);
	}

	FunctionalBasis<1>* CreateLowerDegreeBasis(int degree) override
	{
		return new BernsteinBasis1D(degree);
	}
};

class BernsteinBasis2D : public BernsteinBasis<2>
{
private:
	vector<Bernstein2D> _localFunctions;
private:
	BernsteinBasis2D() {}
public:
	BernsteinBasis2D(int maxPolynomialDegree, bool usePolynomialSpaceQ)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);

		if (usePolynomialSpaceQ)
			Utils::FatalError("Space Q not implemented for the Bernstein basis");
		else
		{
			int nFunctions = Utils::Binomial(maxPolynomialDegree + 2, maxPolynomialDegree);
			_localFunctions.reserve(nFunctions);
			int functionNumber = 0;
			for (int j = 0; j <= maxPolynomialDegree; j++)
			{
				for (int i = 0; i <= maxPolynomialDegree - j; i++)
				{
					_localFunctions.emplace_back(functionNumber, maxPolynomialDegree, i, j);
					functionNumber++;
				}
			}
			assert(nFunctions == _localFunctions.size());
		}
	}

	vector<BasisFunction<2>*> LocalFunctions() override
	{
		vector<BasisFunction<2>*> list;
		list.reserve(_localFunctions.size());
		for (auto& phi : _localFunctions)
			list.push_back(&phi);
		return list;
	}

	int GetDegree() const override
	{
		return _localFunctions[0].GetDegree();
	}

	int Size() const override
	{
		return (int)_localFunctions.size();
	}

	FunctionalBasis<2>* CreateSameBasisForDegree(int degree) override
	{
		return new BernsteinBasis2D(degree, UsePolynomialSpaceQ());
	}

	FunctionalBasis<2>* CreateLowerDegreeBasis(int degree) override
	{
		return new BernsteinBasis2D(degree, UsePolynomialSpaceQ());
	}
};

#ifdef ENABLE_3D

class BernsteinBasis3D : public BernsteinBasis<3>
{
private:
	vector<Bernstein3D> _localFunctions;
private:
	BernsteinBasis3D() {}
public:
	BernsteinBasis3D(int maxPolynomialDegree, bool usePolynomialSpaceQ)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);

		if (usePolynomialSpaceQ)
			Utils::FatalError("Space Q not implemented for the Bernstein basis");
		else
		{
			int nFunctions = Utils::Binomial(maxPolynomialDegree + 3, maxPolynomialDegree);
			_localFunctions.reserve(nFunctions);
			int functionNumber = 0;
			for (int degZ = 0; degZ <= maxPolynomialDegree; degZ++)
			{
				for (int degY = 0; degY <= maxPolynomialDegree - degZ; degY++)
				{
					for (int degX = 0; degX <= maxPolynomialDegree - degY - degZ; degX++)
					{
						_localFunctions.emplace_back(functionNumber, maxPolynomialDegree, degX, degY, degZ);
						functionNumber++;
					}
				}
			}
			assert(nFunctions == _localFunctions.size());
		}
	}

	vector<BasisFunction<3>*> LocalFunctions() override
	{
		vector<BasisFunction<3>*> list;
		list.reserve(_localFunctions.size());
		for (auto& phi : _localFunctions)
			list.push_back(&phi);
		return list;
	}

	int GetDegree() const
	{
		return _localFunctions[0].GetDegree();
	}

	int Size() const override
	{
		return (int)_localFunctions.size();
	}

	FunctionalBasis<3>* CreateSameBasisForDegree(int degree) override
	{
		return new BernsteinBasis3D(degree, UsePolynomialSpaceQ());
	}

	FunctionalBasis<3>* CreateLowerDegreeBasis(int degree) override
	{
		return new BernsteinBasis3D(degree, UsePolynomialSpaceQ());
	}
};
#endif // ENABLE_3D