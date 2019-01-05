#pragma once
#include "Element.h"
#include "ElementInterface.h"

template <class IBasisFunction>
class IPoisson_DGTerms
{
public:
	virtual bool IsGlobalBasis() = 0;

	virtual double VolumicTerm(Element* element, IBasisFunction* func1, IBasisFunction* func2) = 0;

	virtual double CouplingTerm(ElementInterface* interface, Element* element1, IBasisFunction* func1, Element* element2, IBasisFunction* func2) = 0;

	virtual double PenalizationTerm(ElementInterface* interface, Element* element1, IBasisFunction* func1, Element* element2, IBasisFunction* func2, double penalizationCoefficient) = 0;

	virtual double RightHandSide(Element* element, IBasisFunction* func) = 0;
};