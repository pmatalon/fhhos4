#pragma once
#include "../../Utils/Utils.h"


class IDiscreteSpace
{
public:
	// returns [(phi_i|func)]_i
	virtual Vector InnerProdWithBasis(DomFunction func) = 0;

	// returns M^-1 * x
	virtual Vector SolveMassMatrix(const Vector& x) = 0;

	// returns the measure of the support, which equals to (1|1)
	virtual double Measure() = 0;
};
