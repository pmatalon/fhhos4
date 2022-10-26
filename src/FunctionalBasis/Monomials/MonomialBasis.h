#pragma once
#include "../FunctionalBasis.h"
#include "../TensorPolynomial.h"
#include "Monomial1D.h"

template <int Dim>
class MonomialBasis : public FunctionalBasis<Dim>
{
public:
	bool   IsHierarchical()                const override { return true; }
	bool   IsOrthogonalOnCartesianShapes() const override { return false; }
	string BasisCode()                     const override { return "monomials"; }
	static string Code()                                  { return "monomials"; };
};

class MonomialBasis1D : public MonomialBasis<1>
{
private:
	vector<Monomial1D> _localMonomials;
private:
	MonomialBasis1D() {}
public:
	MonomialBasis1D(int maxPolynomialDegree)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);

		_localMonomials.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_localMonomials.emplace_back(i);

	}

	vector<BasisFunction<1>*> LocalFunctions() override
	{
		vector<BasisFunction<1>*> list;
		list.reserve(_localMonomials.size());
		for (Monomial1D& m : _localMonomials)
			list.push_back(&m);
		return list;
	}

	int GetDegree() const override
	{
		return Size() - 1;
	}

	int Size() const override
	{
		return (int)_localMonomials.size();
	}

	bool UsePolynomialSpaceQ() const override
	{
		return false;
	}

	FunctionalBasis<1>* CreateSameBasisForDegree(int degree) override
	{
		return new MonomialBasis1D(degree);
	}

	FunctionalBasis<1>* CreateLowerDegreeBasis(int degree) override
	{
		MonomialBasis1D* lowerBasis = new MonomialBasis1D();
		for (Monomial1D& phi : _localMonomials)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_localMonomials.push_back(phi);
		}
		return lowerBasis;
	}
};

class MonomialBasis2D : public MonomialBasis<2>
{
private:
	vector<Monomial1D> _monomials1D;
	vector<TensorPolynomial2D> _localFunctions;
private:
	MonomialBasis2D() {}
public:
	MonomialBasis2D(int maxPolynomialDegree, bool usePolynomialSpaceQ)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);

		_monomials1D.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_monomials1D.emplace_back(i);

		if (usePolynomialSpaceQ)
		{
			_localFunctions.reserve((maxPolynomialDegree + 1) * (maxPolynomialDegree + 1));
			int functionNumber = 0;
			for (int j = 0; j <= maxPolynomialDegree; j++)
			{
				for (int i = 0; i <= maxPolynomialDegree; i++)
				{
					Monomial1D& polyX = _monomials1D[i];
					Monomial1D& polyY = _monomials1D[j];
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
					Monomial1D& polyX = _monomials1D[i];
					Monomial1D& polyY = _monomials1D[j];
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
		return (int)_monomials1D.size() - 1;
	}

	int Size() const override
	{
		return (int)_localFunctions.size();
	}

	FunctionalBasis<2>* CreateSameBasisForDegree(int degree) override
	{
		return new MonomialBasis2D(degree, UsePolynomialSpaceQ());
	}

	FunctionalBasis<2>* CreateLowerDegreeBasis(int degree) override
	{
		MonomialBasis2D* lowerBasis = new MonomialBasis2D();
		for (Monomial1D& phi : _monomials1D)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_monomials1D.push_back(phi);
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

class MonomialBasis3D : public MonomialBasis<3>
{
private:
	vector<Monomial1D> _monomials1D;
	vector<TensorPolynomial3D> _localFunctions;
private:
	MonomialBasis3D() {}
public:
	MonomialBasis3D(int maxPolynomialDegree, bool usePolynomialSpaceQ)
	{
		maxPolynomialDegree = max(0, maxPolynomialDegree);

		_monomials1D.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_monomials1D.emplace_back(i);

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
						Monomial1D& polyX = _monomials1D[i];
						Monomial1D& polyY = _monomials1D[j];
						Monomial1D& polyZ = _monomials1D[k];
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
						Monomial1D& polyX = _monomials1D[degX];
						Monomial1D& polyY = _monomials1D[degY];
						Monomial1D& polyZ = _monomials1D[degZ];
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

	int GetDegree() const override
	{
		return (int)_monomials1D.size() - 1;
	}

	int Size() const override
	{
		return (int)_localFunctions.size();
	}

	FunctionalBasis<3>* CreateSameBasisForDegree(int degree) override
	{
		return new MonomialBasis3D(degree, UsePolynomialSpaceQ());
	}

	FunctionalBasis<3>* CreateLowerDegreeBasis(int degree) override
	{
		MonomialBasis3D* lowerBasis = new MonomialBasis3D();
		for (Monomial1D& phi : _monomials1D)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_monomials1D.push_back(phi);
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