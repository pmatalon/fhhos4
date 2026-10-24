#pragma once
#include "../FunctionalBasis.h"
#include "../TensorPolynomial.h"
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
	LegendreBasis2D(int maxPolynomialDegree, bool usePolynomialSpaceQ) :
		LegendreBasis2D(usePolynomialSpaceQ ? TensorPartialDegrees2D(max(0, maxPolynomialDegree), 0) : TensorPartialDegrees2D(0, max(0, maxPolynomialDegree)))
	{}

	// The basis of the tensor polynomials of these partial degrees, in their order (see TensorPartialDegrees2D())
	explicit LegendreBasis2D(const vector<array<int, 2>>& partialDegrees)
	{
		int maxPolynomialDegree = 0;
		for (const array<int, 2>& d : partialDegrees)
			maxPolynomialDegree = max({ maxPolynomialDegree, d[0], d[1] });

		_legendre1D.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_legendre1D.emplace_back(i);

		_localFunctions.reserve(partialDegrees.size());
		for (const array<int, 2>& d : partialDegrees)
			_localFunctions.emplace_back((int)_localFunctions.size(), &_legendre1D[d[0]], &_legendre1D[d[1]]);
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

	// The basis of the same space (P or Q) of degree degree
	FunctionalBasis<2>* CreateLowerDegreeBasis(int degree) override
	{
		bool q = UsePolynomialSpaceQ();
		LegendreBasis2D* lowerBasis = new LegendreBasis2D();
		for (Legendre1D& phi : _legendre1D)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_legendre1D.push_back(phi);
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

class LegendreBasis3D : public LegendreBasis<3>
{
private:
	vector<Legendre1D> _legendre1D;
	vector<TensorPolynomial3D> _localFunctions;
private:
	LegendreBasis3D() {}
public:
	LegendreBasis3D(int maxPolynomialDegree, bool usePolynomialSpaceQ) :
		LegendreBasis3D(usePolynomialSpaceQ ? TensorPartialDegrees3D(max(0, maxPolynomialDegree), 0) : TensorPartialDegrees3D(0, max(0, maxPolynomialDegree)))
	{}

	// The basis of the tensor polynomials of these partial degrees, in their order (see TensorPartialDegrees3D())
	explicit LegendreBasis3D(const vector<array<int, 3>>& partialDegrees)
	{
		int maxPolynomialDegree = 0;
		for (const array<int, 3>& d : partialDegrees)
			maxPolynomialDegree = max({ maxPolynomialDegree, d[0], d[1], d[2] });

		_legendre1D.reserve(maxPolynomialDegree + 1);
		for (int i = 0; i <= maxPolynomialDegree; i++)
			_legendre1D.emplace_back(i);

		_localFunctions.reserve(partialDegrees.size());
		for (const array<int, 3>& d : partialDegrees)
			_localFunctions.emplace_back((int)_localFunctions.size(), &_legendre1D[d[0]], &_legendre1D[d[1]], &_legendre1D[d[2]]);
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

	// The basis of the same space (P or Q) of degree degree
	FunctionalBasis<3>* CreateLowerDegreeBasis(int degree) override
	{
		bool q = UsePolynomialSpaceQ();
		LegendreBasis3D* lowerBasis = new LegendreBasis3D();
		for (Legendre1D& phi : _legendre1D)
		{
			if (phi.GetDegree() <= degree)
				lowerBasis->_legendre1D.push_back(phi);
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