#pragma once
#include "Element.h"
#include "IBasisFunction.h"

class IPoisson1D_DGTerms
{
public:
	virtual bool IsGlobalBasis() = 0;

	virtual double VolumicTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2) = 0;

	virtual double MassTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2) = 0;

	virtual double CouplingTerm(BigNumber interface, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2) = 0;

	virtual double PenalizationTerm(BigNumber point, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2, double penalizationCoefficient) = 0;

	virtual double RightHandSide(BigNumber element, IBasisFunction1D* func) = 0;
};