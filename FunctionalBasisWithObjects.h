#pragma once
#include <string>
#include "Element.h"
#include "ElementInterface.h"
//#include "IBasisFunction2D.h"

template <class IBasisFunction>
class FunctionalBasisWithObjects
{
public:
	virtual std::string Name() = 0;

	virtual int GetDegree() = 0;

	virtual int NumberOfLocalFunctionsInElement(Element* element) = 0;

	virtual IBasisFunction* GetLocalBasisFunction(Element* element, int localFunctionNumber) = 0;

	virtual BigNumber GlobalFunctionNumber(Element* element, int localFunctionNumber) = 0;

	virtual double VolumicTerm(Element* element, IBasisFunction* func1, IBasisFunction* func2) = 0;

	virtual double CouplingTerm(ElementInterface* interface, Element* element1, IBasisFunction* func1, Element* element2, IBasisFunction* func2) = 0;

	virtual double PenalizationTerm(ElementInterface* interface, Element* element1, IBasisFunction* func1, Element* element2, IBasisFunction* func2) = 0;

	virtual double RightHandSide(Element* element, IBasisFunction* func) = 0;

	virtual ~FunctionalBasisWithObjects() {}
};