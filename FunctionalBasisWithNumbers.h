#pragma once
#include <string>
#include <map>
#include "Element.h"
#include "IBasisFunction1D.h"
class FunctionalBasisWithNumbers
{
protected:
	std::map<int, IBasisFunction1D*> _localFunctions;

public:
	virtual std::string Name() = 0;

	virtual int GetDegree() = 0;

	int NumberOfLocalFunctionsInElement(BigNumber element)
	{
		return static_cast<int>(this->_localFunctions.size());
	}

	IBasisFunction1D* GetLocalBasisFunction(BigNumber element, int localFunctionNumber)
	{
		return this->_localFunctions[localFunctionNumber];
	}

	BigNumber GlobalFunctionNumber(BigNumber element, int localFunctionNumber)
	{
		return element * NumberOfLocalFunctionsInElement(0) + localFunctionNumber + 1; // +1 so that the numbers start at 1
	}

	virtual double VolumicTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2) = 0;

	virtual double CouplingTerm(BigNumber interface, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2) = 0;

	virtual double PenalizationTerm(BigNumber point, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2, double penalizationCoefficient) = 0;

	virtual double MassTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2) = 0;

	virtual double RightHandSide(BigNumber element, IBasisFunction1D* func) = 0;

	virtual ~FunctionalBasisWithNumbers() {}
};