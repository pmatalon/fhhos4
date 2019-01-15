#pragma once
#include <Eigen/Dense>
#include "IBasisFunction.h"
#include "Element.h"

class Poisson_DG_ReferenceInterface
{
protected:
	Eigen::MatrixXd _couplingTerms;
public:
	Poisson_DG_ReferenceInterface(int nBasisFunctions) :
		_couplingTerms(nBasisFunctions, nBasisFunctions)
	{}

	double CouplingTerm(Element* element1, BasisFunction* phi1, Element* element2, BasisFunction* phi2)
	{
		double result = _couplingTerms(phi1->LocalNumber, phi2->LocalNumber);
		return result;
	}

	virtual void ComputeCouplingTerm(Element* element1, BasisFunction* phi1, Element* element2, BasisFunction* phi2) = 0;
};