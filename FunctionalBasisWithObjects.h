#pragma once
#include <string>
#include "Element.h"
#include "ElementInterface.h"
class FunctionalBasisWithObjects
{
public:
	virtual std::string Name() = 0;

	virtual int NumberOfLocalFunctionsInElement(Element* element) = 0;

	virtual BigNumber GlobalFunctionNumber(Element* element, int localFunctionNumber) = 0;

	virtual double VolumicTerm(Element* element, int localFunctionNumber1, int localFunctionNumber2) = 0;

	virtual double CouplingTerm(ElementInterface* interface, Element* element1, int localFunctionNumber1, Element* element2, int localFunctionNumber2) = 0;

	virtual double PenalizationTerm(ElementInterface* interface, Element* element1, int localFunctionNumber1, Element* element2, int localFunctionNumber2) = 0;

	virtual double RightHandSide(Element* element, int localFunctionNumber) = 0;

	virtual ~FunctionalBasisWithObjects() {}
};