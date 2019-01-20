#pragma once
#include "Element.h"
#include "IBasisFunction.h"

class Poisson_DG_Face
{
public:
	virtual double CouplingTerm(Poisson_DG_Element* element1, BasisFunction* phi1, Poisson_DG_Element* element2, BasisFunction* phi2) = 0;
	virtual double PenalizationTerm(Poisson_DG_Element* element1, BasisFunction* phi1, Poisson_DG_Element* element2, BasisFunction* phi2, double penalizationCoefficient) = 0;
};