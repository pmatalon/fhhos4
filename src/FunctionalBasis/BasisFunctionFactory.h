#pragma once
#include <cstdio>
#include "Monomial1D.h"
#include "Legendre1D.h"
#include "Bernstein1D.h"
#include "Hemker1D.h"
using namespace std;

class BasisFunctionFactory
{
public:
	static IBasisFunction1D* Create(string basisCode, int maxPolynomialDegree, int i)
	{
		if (basisCode.compare(Monomial1D::Code()) == 0)
			return new Monomial1D(i);
		if (basisCode.compare(Legendre1D::Code()) == 0)
			return new Legendre1D(i);
		if (basisCode.compare(NormalizedLegendre1D::Code()) == 0)
			return new NormalizedLegendre1D(i);
		if (basisCode.compare(Bernstein1D::Code()) == 0)
			return new Bernstein1D(maxPolynomialDegree, i);
		if (basisCode.compare(Hemker1D::Code()) == 0)
			return new Hemker1D(maxPolynomialDegree, i);

		cout << "Basis '" << basisCode << "' not managed!";
		exit(EXIT_FAILURE);
		return NULL;
	}

	static bool IsHierarchicalBasis(string basisCode)
	{
		return basisCode.compare(Monomial1D::Code()) == 0 || basisCode.compare(Legendre1D::Code()) == 0;
	}

	static bool IsOrthogonalBasis(string basisCode)
	{
		return basisCode.compare(Legendre1D::Code()) == 0;
	}
};
