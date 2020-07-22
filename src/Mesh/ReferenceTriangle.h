#pragma once
#include "ReferenceShape.h"
#include "../Utils/Dunavant/Dunavant.h"

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

	inline double Diameter() const override
	{
		assert(false);
	}
	inline double Measure() const override
	{
		return _measure;
	}
	inline DomPoint Center() const override
	{
		assert(false);
	}

	vector<RefPoint> QuadraturePoints() const
	{
		Dunavant dunavant;
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