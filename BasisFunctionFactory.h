#pragma once
#include "MonomialBasis1D.h"
#include "LegendreBasis1D.h"
#include "BernsteinBasis1D.h"
#include "Bernstein2Basis1D.h"

class BasisFunctionFactory
{
public:
	static IBasisFunction1D* Create(std::string polynomialCode, int maxPolynomialDegree, int i)
	{
		if (polynomialCode.compare(Monomial1D::Code()) == 0)
			return new Monomial1D(i);
		if (polynomialCode.compare(Legendre1D::Code()) == 0)
			return new Legendre1D(i, false);
		if (polynomialCode.compare(Bernstein1D::Code()) == 0)
			return new Bernstein1D(maxPolynomialDegree, i);
		if (polynomialCode.compare(Bernstein2_1D::Code()) == 0)
			return new Bernstein2_1D(maxPolynomialDegree, i);
		return NULL;
	}
};