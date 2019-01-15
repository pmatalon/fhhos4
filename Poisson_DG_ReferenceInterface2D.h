#pragma once
#include "Poisson_DG_ReferenceInterface.h"

class Poisson_DG_ReferenceInterface2D : public Poisson_DG_ReferenceInterface
{
public:
	Poisson_DG_ReferenceInterface2D(int nBasisFunctions) :
		Poisson_DG_ReferenceInterface(nBasisFunctions)
	{}

	void ComputeCouplingTerm(Element* element1, BasisFunction* p_phi1, Element* element2, BasisFunction* p_phi2)
	{
		IBasisFunction2D* phi1 = dynamic_cast<IBasisFunction2D*>(p_phi1);
		IBasisFunction2D* phi2 = dynamic_cast<IBasisFunction2D*>(p_phi2);


	}
};