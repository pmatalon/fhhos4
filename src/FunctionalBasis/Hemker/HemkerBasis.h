#pragma once
#include "../FunctionalBasis.h"
#include "../TensorPolynomial.h"
#include "Hemker1D.h"

template <int Dim>
class HemkerBasis : public FunctionalBasis<Dim>
{
public:
	bool   IsHierarchical()                const override { return false; }
	string BasisCode()                     const override { return "hemker"; }
	static string Code()                                  { return "hemker"; };
};

class HemkerBasis1D : public HemkerBasis<1>
{
private:
	vector<Hemker1D> _localFunctions;
private:
	HemkerBasis1D() {}
public:
	HemkerBasis1D(int maxPolynomialDegree)
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
		for (Hemker1D& m : _localFunctions)
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
		return new HemkerBasis1D(degree);
	}

	FunctionalBasis<1>* CreateLowerDegreeBasis(int degree) override
	{
		return new HemkerBasis1D(degree);
	}
};

class HemkerBasis2D : public HemkerBasis<2>
{
private:
	vector<Hemker1D> _hemker1D;
	vector<TensorPolynomial2D> _localFunctions;
private:
	HemkerBasis2D() {}
public:
	HemkerBasis2D(int maxPolynomialDegree, bool usePolynomialSpaceQ)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);

		_hemker1D.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_hemker1D.emplace_back(maxPolynomialDegree, i);

		if (usePolynomialSpaceQ)
		{
			_localFunctions.reserve((maxPolynomialDegree + 1) * (maxPolynomialDegree + 1));
			int functionNumber = 0;
			for (int j = 0; j <= maxPolynomialDegree; j++)
			{
				for (int i = 0; i <= maxPolynomialDegree; i++)
				{
					Hemker1D& polyX = _hemker1D[i];
					Hemker1D& polyY = _hemker1D[j];
					_localFunctions.emplace_back(functionNumber, &polyX, &polyY);
					functionNumber++;
				}
			}
		}
		else
		{
			int nFunctions = Utils::Binomial(maxPolynomialDegree + 2, maxPolynomialDegree);
			_localFunctions.reserve(nFunctions);
			int functionNumber = 0;
			for (int degree = 0; degree <= maxPolynomialDegree; degree++)
			{
				for (int j = 0; j <= degree; j++)
				{
					int i = degree - j;
					Hemker1D* polyX = new Hemker1D(i, i); // TODO: free
					Hemker1D* polyY = new Hemker1D(j, j); // TODO: free
					_localFunctions.emplace_back(functionNumber, polyX, polyY);
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
		for (TensorPolynomial2D& tp : _localFunctions)
			list.push_back(&tp);
		return list;
	}

	int GetDegree() const override
	{
		return (int)_hemker1D.size() - 1;
	}

	int Size() const override
	{
		return (int)_localFunctions.size();
	}

	FunctionalBasis<2>* CreateSameBasisForDegree(int degree) override
	{
		return new HemkerBasis2D(degree, UsePolynomialSpaceQ());
	}

	FunctionalBasis<2>* CreateLowerDegreeBasis(int degree) override
	{
		return new HemkerBasis2D(degree, UsePolynomialSpaceQ());
	}
};

#ifdef ENABLE_3D

class HemkerBasis3D : public HemkerBasis<3>
{
private:
	vector<Hemker1D> _hemker1D;
	vector<TensorPolynomial3D> _localFunctions;
private:
	HemkerBasis3D() {}
public:
	HemkerBasis3D(int maxPolynomialDegree, bool usePolynomialSpaceQ)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);

		_hemker1D.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_hemker1D.emplace_back(maxPolynomialDegree, i);

		if (usePolynomialSpaceQ)
		{
			_localFunctions.reserve(pow(maxPolynomialDegree + 1, 3));
			int functionNumber = 0;
			for (int k = 0; k <= maxPolynomialDegree; k++)
			{
				for (int j = 0; j <= maxPolynomialDegree; j++)
				{
					for (int i = 0; i <= maxPolynomialDegree; i++)
					{
						Hemker1D& polyX = _hemker1D[i];
						Hemker1D& polyY = _hemker1D[j];
						Hemker1D& polyZ = _hemker1D[k];
						_localFunctions.emplace_back(functionNumber, &polyX, &polyY, &polyZ);
						functionNumber++;
					}
				}
			}
		}
		else
		{
			int nFunctions = Utils::Binomial(maxPolynomialDegree + 3, maxPolynomialDegree);
			_localFunctions.reserve(nFunctions);
			int functionNumber = 0;
			for (int degree = 0; degree <= maxPolynomialDegree; degree++)
			{
				for (int k = 0; k <= degree; k++)
				{
					for (int j = 0; j <= degree - k; j++)
					{
						int i = degree - k - j;
						Hemker1D* polyX = new Hemker1D(i, i); // TODO: free
						Hemker1D* polyY = new Hemker1D(j, j); // TODO: free
						Hemker1D* polyZ = new Hemker1D(k, k); // TODO: free
						_localFunctions.emplace_back(functionNumber, polyX, polyY, polyZ);
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
		for (TensorPolynomial3D& tp : _localFunctions)
			list.push_back(&tp);
		return list;
	}

	int GetDegree() const
	{
		return (int)_hemker1D.size() - 1;
	}

	int Size() const override
	{
		return (int)_localFunctions.size();
	}

	FunctionalBasis<3>* CreateSameBasisForDegree(int degree) override
	{
		return new HemkerBasis3D(degree, UsePolynomialSpaceQ());
	}

	FunctionalBasis<3>* CreateLowerDegreeBasis(int degree) override
	{
		return new HemkerBasis3D(degree, UsePolynomialSpaceQ());
	}
};
#endif // ENABLE_3D