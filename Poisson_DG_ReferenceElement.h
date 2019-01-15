#pragma once
#include <Eigen/Dense>
#include "IBasisFunction.h"

class Poisson_DG_ReferenceElement
{
protected:
	Eigen::MatrixXd _volumicTerms;
public:
	Poisson_DG_ReferenceElement(int nBasisFunctions) :
		_volumicTerms(nBasisFunctions, nBasisFunctions)
	{}

	double VolumicTerm(BasisFunction* phi1, BasisFunction* phi2)
	{
		double result = _volumicTerms(phi1->LocalNumber, phi2->LocalNumber);
		return result;
	}

	virtual void ComputeVolumicTerm(BasisFunction* phi1, BasisFunction* phi2) = 0;
};