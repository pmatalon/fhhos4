#pragma once
#include "Element.h"
#include "IBasisFunction.h"

class Poisson_DG_ReferenceFace;

class Poisson_DG_Face
{
public:
	virtual double CouplingTerm(Element* element1, BasisFunction* phi1, Element* element2, BasisFunction* phi2, Poisson_DG_ReferenceFace* referenceInterface) = 0;
};