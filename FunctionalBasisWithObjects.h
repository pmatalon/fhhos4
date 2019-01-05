#pragma once
#include <string>
#include <map>
#include "Element.h"
#include "ElementInterface.h"

template <class IBasisFunction>
class FunctionalBasisWithObjects
{
protected:
	map<int, IBasisFunction*> _localFunctions;

public:
	virtual std::string Name() = 0;

	virtual int GetDegree() = 0;

	int NumberOfLocalFunctionsInElement(Element* element)
	{
		return static_cast<int>(this->_localFunctions.size());
	}

	IBasisFunction* GetLocalBasisFunction(Element* element, int localFunctionNumber)
	{
		return this->_localFunctions[localFunctionNumber];
	}

	BigNumber GlobalFunctionNumber(Element* element, int localFunctionNumber)
	{
		return element->Number * static_cast<int>(this->_localFunctions.size()) + localFunctionNumber + 1; // +1 so that the numbers start at 1
	}

	virtual double VolumicTerm(Element* element, IBasisFunction* func1, IBasisFunction* func2) = 0;

	virtual double CouplingTerm(ElementInterface* interface, Element* element1, IBasisFunction* func1, Element* element2, IBasisFunction* func2) = 0;

	virtual double PenalizationTerm(ElementInterface* interface, Element* element1, IBasisFunction* func1, Element* element2, IBasisFunction* func2, double penalizationCoefficient) = 0;

	virtual double RightHandSide(Element* element, IBasisFunction* func) = 0;

	virtual ~FunctionalBasisWithObjects() {}
};