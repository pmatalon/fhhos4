/*#pragma once
#include "Poisson_DG_Element.h"
#include "Square.h"
#include "IBasisFunction.h"

class Poisson_DG_Square : public Poisson_DG_Element
{
private:
	Square* _square;
public:
	Poisson_DG_Square(Square* square)
	{
		this->_square = square;
	}

	double VolumicTerm(BasisFunction* phi1, BasisFunction* phi2, Poisson_DG_ReferenceElement* referenceElement)
	{
		return referenceElement->VolumicTerm(phi1, phi2);
	}
};*/