#pragma once
#include "Poisson_DG_Interface.h"
#include "Element2DInterface.h"
#include "IBasisFunction.h"
#include "Poisson_DG_ReferenceInterface.h"

class Poisson_DG_Interface2D : public Poisson_DG_Interface
{
private:
	Element2DInterface* _interface;
public:
	Poisson_DG_Interface2D(Element2DInterface* interface)
	{
		this->_interface = interface;
	}

	double VolumicTerm(Element* element1, BasisFunction* phi1, Element* element2, BasisFunction* phi2, Poisson_DG_ReferenceInterface* referenceInterface)
	{
		return referenceInterface->CouplingTerm(element1, phi1, element2, phi2);
	}
};