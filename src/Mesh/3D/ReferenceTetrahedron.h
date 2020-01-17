#pragma once
#include "../ReferenceShape.h"
#include "../../Utils/Keast/Keast.h"

class ReferenceTetrahedron : public ReferenceShape<3>
{
private:
	RefPoint A;
	RefPoint B;
	RefPoint C;
	RefPoint D;
	double _measure;

public:
	ReferenceTetrahedron() :
		ReferenceShape<3>(),
		A(0, 0, 0), B(1, 0, 0), C(0, 1, 0), D(0, 0, 1)
	{
		DimMatrix<3> m;
		m.col(0) = Vect<3>(A, B);
		m.col(1) = Vect<3>(A, C);
		m.col(2) = Vect<3>(A, D);
		_measure = abs(m.determinant()) / 6;
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

	double Integral(RefFunction func) const override
	{
		Keast keast;
		return _measure * keast.Quadrature(func);
	}
	double Integral(RefFunction func, int polynomialDegree) const override
	{
		Keast keast(polynomialDegree);
		return _measure * keast.Quadrature(func);
	}

};