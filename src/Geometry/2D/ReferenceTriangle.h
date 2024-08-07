#pragma once
#include "../ReferenceShape.h"
#include "../../QuadratureRules/Dunavant/Dunavant.h"

class ReferenceTriangle : public ReferenceShape<2>
{
private:
	RefPoint A;
	RefPoint B;
	RefPoint C;
	double _measure;

public:
	ReferenceTriangle() : 
		ReferenceShape<2>(),
		A(0, 0), B(1, 0), C(0, 1)
	{
		_measure = 0.5 * abs(A.X * (B.Y - C.Y) + B.X * (C.Y - A.Y) + C.X * (A.Y - B.Y));
	}

	string Name() const override { return "Reference Triangle"; }

	inline double Diameter() const override
	{
		assert(false);
		return 0.0; // to avoid warning
	}
	inline double Measure() const override
	{
		return _measure;
	}
	inline DomPoint Center() const override
	{
		assert(false);
		return DomPoint(); // to avoid warning
	}

	vector<RefPoint> QuadraturePoints() const override
	{
		Dunavant dunavant;
		return dunavant.Points();
	}

	vector<RefPoint> QuadraturePoints(int polynomialDegree) const override
	{
		Dunavant dunavant(polynomialDegree);
		return dunavant.Points();
	}

	double Integral(RefFunction func) const override
	{
		Dunavant dunavant;
		return _measure * dunavant.Quadrature(func);
	}
	double Integral(RefFunction func, int polynomialDegree) const override
	{
		Dunavant dunavant(polynomialDegree);
		return _measure * dunavant.Quadrature(func);
	}

};