#pragma once
#include "FunctionalBasis.h"
#include "TensorPolynomial.h"
#include "Legendre1D.h"

template <int Dim>
class LegendreBasis : public FunctionalBasis<Dim>
{
public:
	bool   IsHierarchical()                const override { return true; }
	bool   IsOrthogonalOnCartesianShapes() const override { return true; }
	string BasisCode()                     const override { return "legendre"; }
	static string Code()                                  { return "legendre"; };
};

class LegendreBasis1D : public LegendreBasis<1>
{
private:
	vector<Legendre1D> _localFunctions;
private:
	LegendreBasis1D() {}
public:
	LegendreBasis1D(int maxPolynomialDegree)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);

		_localFunctions.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_localFunctions.emplace_back(i);

	}
	vector<BasisFunction<1>*> LocalFunctions() override
	{
		vector<BasisFunction<1>*> list;
		list.reserve(_localFunctions.size());
		for (Legendre1D& m : _localFunctions)
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
		return new LegendreBasis1D(degree);
	}

	FunctionalBasis<1>* CreateLowerDegreeBasis(int degree) override
	{
		LegendreBasis1D* lowerBasis = new LegendreBasis1D();
		for (Legendre1D& phi : _localFunctions)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_localFunctions.push_back(phi);
		}
		return lowerBasis;
	}
};

class LegendreBasis2D : public LegendreBasis<2>
{
private:
	vector<Legendre1D> _legendre1D;
	vector<TensorPolynomial2D> _localFunctions;
private:
	LegendreBasis2D() {}
public:
	LegendreBasis2D(int maxPolynomialDegree, bool usePolynomialSpaceQ)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);

		_legendre1D.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_legendre1D.emplace_back(i);

		if (usePolynomialSpaceQ)
		{
			_localFunctions.reserve((maxPolynomialDegree + 1) * (maxPolynomialDegree + 1));
			int functionNumber = 0;
			for (int j = 0; j <= maxPolynomialDegree; j++)
			{
				for (int i = 0; i <= maxPolynomialDegree; i++)
				{
					Legendre1D& polyX = _legendre1D[i];
					Legendre1D& polyY = _legendre1D[j];
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
					Legendre1D& polyX = _legendre1D[i];
					Legendre1D& polyY = _legendre1D[j];
					_localFunctions.emplace_back(functionNumber, &polyX, &polyY);
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
		return (int)_legendre1D.size() - 1;
	}

	int Size() const override
	{
		return (int)_localFunctions.size();
	}

	FunctionalBasis<2>* CreateSameBasisForDegree(int degree) override
	{
		return new LegendreBasis2D(degree, UsePolynomialSpaceQ());
	}

	FunctionalBasis<2>* CreateLowerDegreeBasis(int degree) override
	{
		LegendreBasis2D* lowerBasis = new LegendreBasis2D();
		for (Legendre1D& phi : _legendre1D)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_legendre1D.push_back(phi);
		}
		for (TensorPolynomial2D& tp : _localFunctions)
		{
			if (tp.GetDegree() <= degree)
				lowerBasis->_localFunctions.push_back(tp);
		}
		return lowerBasis;
	}
};

#ifdef ENABLE_3D

class LegendreBasis3D : public LegendreBasis<3>
{
private:
	vector<Legendre1D> _legendre1D;
	vector<TensorPolynomial3D> _localFunctions;
private:
	LegendreBasis3D() {}
public:
	LegendreBasis3D(int maxPolynomialDegree, bool usePolynomialSpaceQ)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);

		_legendre1D.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_legendre1D.emplace_back(i);

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
						Legendre1D& polyX = _legendre1D[i];
						Legendre1D& polyY = _legendre1D[j];
						Legendre1D& polyZ = _legendre1D[k];
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
				for (int degZ = 0; degZ <= degree; degZ++)
				{
					for (int degY = 0; degY <= degree - degZ; degY++)
					{
						int degX = degree - degZ - degY;
						Legendre1D& polyX = _legendre1D[degX];
						Legendre1D& polyY = _legendre1D[degY];
						Legendre1D& polyZ = _legendre1D[degZ];
						_localFunctions.emplace_back(functionNumber, &polyX, &polyY, &polyZ);
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
		return (int)_legendre1D.size() - 1;
	}

	int Size() const override
	{
		return (int)_localFunctions.size();
	}

	FunctionalBasis<3>* CreateSameBasisForDegree(int degree) override
	{
		return new LegendreBasis3D(degree, UsePolynomialSpaceQ());
	}

	FunctionalBasis<3>* CreateLowerDegreeBasis(int degree) override
	{
		LegendreBasis3D* lowerBasis = new LegendreBasis3D();
		for (Legendre1D& phi : _legendre1D)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_legendre1D.push_back(phi);
		}
		for (TensorPolynomial3D& tp : _localFunctions)
		{
			if (tp.GetDegree() <= degree)
				lowerBasis->_localFunctions.push_back(tp);
		}
		return lowerBasis;
	}
};
#endif // ENABLE_3D