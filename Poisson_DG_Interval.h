/*#pragma once
#include "Poisson_DG_Element.h"
#include "Interval.h"
#include "IBasisFunction.h"

class Poisson_DG_Interval : public Poisson_DG_Element
{
private:
	Interval* _interval;
public:
	Poisson_DG_Interval(Interval* interval)
	{
		this->_interval = interval;
	}

	double VolumicTerm(BasisFunction* phi1, BasisFunction* phi2, Poisson_DG_ReferenceElement* referenceElement)
	{
		double h = this->_interval->B - this->_interval->A;
		double factor = 2 / h;
		return factor * referenceElement->VolumicTerm(phi1, phi2);
	}
};*/