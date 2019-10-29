#pragma once
#include "ReferenceElement.h"
#include "../Utils/Dunavant/Dunavant.h"

class ReferenceTriangle : public ReferenceElement<2>
{
private:
	RefPoint A;
	RefPoint B;
	RefPoint C;
	double _measure;

public:
	ReferenceTriangle() : 
		ReferenceElement<2>(),
		A(0, 0), B(1, 0), C(0, 1)
	{
		_measure = 0.5 * abs(A.X * (B.Y - C.Y) + B.X * (C.Y - A.Y) + C.X * (A.Y - B.Y));
	}

	inline double Measure()
	{
		return _measure;
	}

	double ComputeIntegral(RefFunction func) const override
	{
		Dunavant dunavant;
		return _measure * dunavant.Quadrature(func);
	}
	double ComputeIntegral(RefFunction func, int polynomialDegree) const override
	{
		Dunavant dunavant(polynomialDegree);
		return _measure * dunavant.Quadrature(func);
	}
	double ComputeIntegral(BasisFunction<2>* phi) const
	{
		return ReferenceElement<2>::ComputeIntegral(phi);
	}

};