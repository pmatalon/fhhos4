#pragma once
#include <string>
#include <map>
#include "Element.h"
#include "IBasisFunction.h"
class FunctionalBasisWithNumbers
{
protected:
	std::map<int, IBasisFunction1D*> LocalFunctions;

public:
	virtual std::string Name() = 0;

	virtual int GetDegree() = 0;

	int NumberOfLocalFunctionsInElement(BigNumber element)
	{
		return static_cast<int>(this->LocalFunctions.size());
	}

	IBasisFunction1D* GetLocalBasisFunction(BigNumber element, int localFunctionNumber)
	{
		return this->LocalFunctions[localFunctionNumber];
	}

	BigNumber GlobalFunctionNumber(BigNumber element, int localFunctionNumber)
	{
		return element * NumberOfLocalFunctionsInElement(0) + localFunctionNumber; // the numbers start at 0
	}

	virtual ~FunctionalBasisWithNumbers() {}
};