#pragma once
#include <cstdio>
#include "Monomials/MonomialBasis.h"
#include "Legendre/LegendreBasis.h"
#include "Bernstein/BernsteinBasis.h"
#include "Hemker/HemkerBasis.h"
#include "Lagrange/LagrangeBasis.h"
using namespace std;

template <int Dim>
class FunctionalBasisFactory
{
public:
	static FunctionalBasis<Dim>* Create(string basisCode, int maxPolynomialDegree, bool usePolynomialSpaceQ) { assert(false); }
	// The basis of Q^qDegree + P^pDegree (see TensorPolynomial.h), in 2D and 3D: monomials and Legendre only
	static FunctionalBasis<Dim>* CreateQPlusP(string basisCode, int qDegree, int pDegree) { assert(false); }
};


template <>
inline FunctionalBasis<0>* FunctionalBasisFactory<0>::Create(string basisCode, int maxPolynomialDegree, bool usePolynomialSpaceQ)
{
	return new FunctionalBasis0D();
}

template <>
inline FunctionalBasis<1>* FunctionalBasisFactory<1>::Create(string basisCode, int maxPolynomialDegree, bool usePolynomialSpaceQ)
{
	if (basisCode.compare(MonomialBasis<1>::Code()) == 0)
		return new MonomialBasis1D(maxPolynomialDegree);
	if (basisCode.compare(LegendreBasis<1>::Code()) == 0)
		return new LegendreBasis1D(maxPolynomialDegree); 
	if (basisCode.compare(HemkerBasis<1>::Code()) == 0)
		return new HemkerBasis1D(maxPolynomialDegree);
	if (basisCode.compare(BernsteinBasis<1>::Code()) == 0)
		return new BernsteinBasis1D(maxPolynomialDegree);

	Utils::FatalError("Basis '" + basisCode + "' not managed!");
	return nullptr;
}

template <>
inline FunctionalBasis<1>* FunctionalBasisFactory<1>::CreateQPlusP(string basisCode, int qDegree, int pDegree)
{
	return Create(basisCode, max(qDegree, pDegree), false); // Q^d = P^d in 1D
}

template <>
inline FunctionalBasis<2>* FunctionalBasisFactory<2>::Create(string basisCode, int maxPolynomialDegree, bool usePolynomialSpaceQ)
{
	if (basisCode.compare(MonomialBasis<2>::Code()) == 0)
		return new MonomialBasis2D(maxPolynomialDegree, usePolynomialSpaceQ);
	if (basisCode.compare(LegendreBasis<2>::Code()) == 0)
		return new LegendreBasis2D(maxPolynomialDegree, usePolynomialSpaceQ);
	if (basisCode.compare(HemkerBasis<2>::Code()) == 0)
		return new HemkerBasis2D(maxPolynomialDegree, usePolynomialSpaceQ);
	if (basisCode.compare(BernsteinBasis<2>::Code()) == 0)
		return new BernsteinBasis2D(maxPolynomialDegree, usePolynomialSpaceQ);
	if (basisCode.compare(LagrangeBasis<2>::Code()) == 0)
		return new LagrangeP1Basis(maxPolynomialDegree);

	Utils::FatalError("Basis '" + basisCode + "' not managed!");
	return nullptr;
}

template <>
inline FunctionalBasis<2>* FunctionalBasisFactory<2>::CreateQPlusP(string basisCode, int qDegree, int pDegree)
{
	if (basisCode.compare(MonomialBasis<2>::Code()) == 0)
		return new MonomialBasis2D(TensorPartialDegrees2D(max(0, qDegree), pDegree));
	if (basisCode.compare(LegendreBasis<2>::Code()) == 0)
		return new LegendreBasis2D(TensorPartialDegrees2D(max(0, qDegree), pDegree));

	Utils::FatalError("Basis '" + basisCode + "': the space Q + P is only implemented with the bases monomials and legendre.");
	return nullptr;
}

#ifdef ENABLE_3D

template <>
inline FunctionalBasis<3>* FunctionalBasisFactory<3>::Create(string basisCode, int maxPolynomialDegree, bool usePolynomialSpaceQ)
{
	if (basisCode.compare(MonomialBasis<3>::Code()) == 0)
		return new MonomialBasis3D(maxPolynomialDegree, usePolynomialSpaceQ);
	if (basisCode.compare(LegendreBasis<3>::Code()) == 0)
		return new LegendreBasis3D(maxPolynomialDegree, usePolynomialSpaceQ);
	if (basisCode.compare(HemkerBasis<3>::Code()) == 0)
		return new HemkerBasis3D(maxPolynomialDegree, usePolynomialSpaceQ);
	if (basisCode.compare(BernsteinBasis<3>::Code()) == 0)
		return new BernsteinBasis3D(maxPolynomialDegree, usePolynomialSpaceQ);

	Utils::FatalError("Basis '" + basisCode + "' not managed!");
	return nullptr;
}

template <>
inline FunctionalBasis<3>* FunctionalBasisFactory<3>::CreateQPlusP(string basisCode, int qDegree, int pDegree)
{
	if (basisCode.compare(MonomialBasis<3>::Code()) == 0)
		return new MonomialBasis3D(TensorPartialDegrees3D(max(0, qDegree), pDegree));
	if (basisCode.compare(LegendreBasis<3>::Code()) == 0)
		return new LegendreBasis3D(TensorPartialDegrees3D(max(0, qDegree), pDegree));

	Utils::FatalError("Basis '" + basisCode + "': the space Q + P is only implemented with the bases monomials and legendre.");
	return nullptr;
}

#endif // ENABLE_3D
