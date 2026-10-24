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
	MonomialBasis2D(int maxPolynomialDegree, bool usePolynomialSpaceQ) :
		MonomialBasis2D(usePolynomialSpaceQ ? TensorPartialDegrees2D(max(0, maxPolynomialDegree), 0) : TensorPartialDegrees2D(0, max(0, maxPolynomialDegree)))
	{}

	// The basis of the tensor polynomials of these partial degrees, in their order (see TensorPartialDegrees2D())
	explicit MonomialBasis2D(const vector<array<int, 2>>& partialDegrees)
	{
		int maxPolynomialDegree = 0;
		for (const array<int, 2>& d : partialDegrees)
			maxPolynomialDegree = max({ maxPolynomialDegree, d[0], d[1] });

		_monomials1D.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_monomials1D.emplace_back(i);

		_localFunctions.reserve(partialDegrees.size());
		for (const array<int, 2>& d : partialDegrees)
			_localFunctions.emplace_back((int)_localFunctions.size(), &_monomials1D[d[0]], &_monomials1D[d[1]]);
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

	// The basis of the same space (P or Q) of degree degree
	FunctionalBasis<2>* CreateLowerDegreeBasis(int degree) override
	{
		bool q = UsePolynomialSpaceQ();
		MonomialBasis2D* lowerBasis = new MonomialBasis2D();
		for (Monomial1D& phi : _monomials1D)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_monomials1D.push_back(phi);
		}
		for (TensorPolynomial2D& tp : _localFunctions)
		{
			if ((q ? tp.GetQDegree() : tp.GetDegree()) <= degree)
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
	MonomialBasis3D(int maxPolynomialDegree, bool usePolynomialSpaceQ) :
		MonomialBasis3D(usePolynomialSpaceQ ? TensorPartialDegrees3D(max(0, maxPolynomialDegree), 0) : TensorPartialDegrees3D(0, max(0, maxPolynomialDegree)))
	{}

	// The basis of the tensor polynomials of these partial degrees, in their order (see TensorPartialDegrees3D())
	explicit MonomialBasis3D(const vector<array<int, 3>>& partialDegrees)
	{
		int maxPolynomialDegree = 0;
		for (const array<int, 3>& d : partialDegrees)
			maxPolynomialDegree = max({ maxPolynomialDegree, d[0], d[1], d[2] });

		_monomials1D.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_monomials1D.emplace_back(i);

		_localFunctions.reserve(partialDegrees.size());
		for (const array<int, 3>& d : partialDegrees)
			_localFunctions.emplace_back((int)_localFunctions.size(), &_monomials1D[d[0]], &_monomials1D[d[1]], &_monomials1D[d[2]]);
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

	// The basis of the same space (P or Q) of degree degree
	FunctionalBasis<3>* CreateLowerDegreeBasis(int degree) override
	{
		bool q = UsePolynomialSpaceQ();
		MonomialBasis3D* lowerBasis = new MonomialBasis3D();
		for (Monomial1D& phi : _monomials1D)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_monomials1D.push_back(phi);
		}
		for (TensorPolynomial3D& tp : _localFunctions)
		{
			if ((q ? tp.GetQDegree() : tp.GetDegree()) <= degree)
				lowerBasis->_localFunctions.push_back(tp);
		}
		return lowerBasis;
	}
};
#endif // ENABLE_3D