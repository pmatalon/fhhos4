#pragma once
#include <string>
#include "Element.h"
class FunctionalBasisWithNumbers
{
public:
	virtual std::string Name() = 0;

	virtual int NumberOfLocalFunctionsInElement(BigNumber element) = 0;

	virtual BigNumber GlobalFunctionNumber(BigNumber element, int localFunctionNumber) = 0;

	virtual double VolumicTerm(BigNumber element, int localFunctionNumber1, int localFunctionNumber2) = 0;

	virtual double CouplingTerm(BigNumber interface, BigNumber element1, int localFunctionNumber1, BigNumber element2, int localFunctionNumber2) = 0;

	virtual double PenalizationTerm(BigNumber point, BigNumber element1, int localFunctionNumber1, BigNumber element2, int localFunctionNumber2) = 0;

	virtual double RightHandSide(BigNumber element, int localFunctionNumber) = 0;

	virtual ~FunctionalBasisWithNumbers() {}
};