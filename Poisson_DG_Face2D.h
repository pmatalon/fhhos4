/*#pragma once
#include "Poisson_DG_Face.h"
#include "Face2D.h"
#include "IBasisFunction.h"
#include "Poisson_DG_ReferenceFace.h"

class Poisson_DG_Face2D : public Poisson_DG_Face
{
private:
	Face2D* _interface;
public:
	Poisson_DG_Face2D(Face2D* interface)
	{
		this->_interface = interface;
	}

	double CouplingTerm(Element* element1, BasisFunction* phi1, Element* element2, BasisFunction* phi2, Poisson_DG_ReferenceFace* referenceInterface)
	{
		return referenceInterface->CouplingTerm(element1, phi1, element2, phi2);
	}
};*/