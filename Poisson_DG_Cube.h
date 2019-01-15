#pragma once
#include "Poisson_DG_Element.h"
#include "Cube.h"
#include "IBasisFunction.h"

class Poisson_DG_Cube : public Poisson_DG_Element
{
private:
	Cube* _cube;
public:
	Poisson_DG_Cube(Cube* cube)
	{
		this->_cube = cube;
	}

	double VolumicTerm(BasisFunction* phi1, BasisFunction* phi2, Poisson_DG_ReferenceElement* referenceElement)
	{
		double h = this->_cube->Width;
		DefInterval refInterval = phi1->DefinitionInterval();
		double factor = h / refInterval.Length;
		return factor * referenceElement->VolumicTerm(phi1, phi2);
	}
};