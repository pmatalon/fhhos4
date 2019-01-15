#pragma once
#include "Element.h"
#include "IBasisFunction.h"

class Poisson_DG_ReferenceInterface;

class Poisson_DG_Interface
{
public:
	virtual double CouplingTerm(Element* element1, BasisFunction* phi1, Element* element2, BasisFunction* phi2, Poisson_DG_ReferenceInterface* referenceInterface) = 0;
};