#pragma once
#include <string>
#include "Element.h"
#include "IBasisFunction1D.h"
class FunctionalBasisWithNumbers
{
public:
	virtual std::string Name() = 0;

	virtual int NumberOfLocalFunctionsInElement(BigNumber element) = 0;

	virtual IBasisFunction1D* GetLocalBasisFunction(BigNumber element, int localFunctionNumber) = 0;

	virtual BigNumber GlobalFunctionNumber(BigNumber element, int localFunctionNumber) = 0;

	virtual double VolumicTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2) = 0;

	virtual double CouplingTerm(BigNumber interface, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2) = 0;

	virtual double PenalizationTerm(BigNumber point, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2) = 0;

	virtual double RightHandSide(BigNumber element, IBasisFunction1D* func) = 0;

	virtual ~FunctionalBasisWithNumbers() {}
};