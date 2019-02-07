#pragma once
#include <Eigen/Dense>
#include "BasisFunction.h"

class Poisson_DG_ReferenceElement
{
protected:
	Eigen::MatrixXd _volumicTerms;
	Eigen::MatrixXd _massTerms;
public:
	Poisson_DG_ReferenceElement(int nBasisFunctions) :
		_volumicTerms(nBasisFunctions, nBasisFunctions), _massTerms(nBasisFunctions, nBasisFunctions)
	{}

	double VolumicTerm(BasisFunction* phi1, BasisFunction* phi2)
	{
		return this->_volumicTerms(phi1->LocalNumber, phi2->LocalNumber);
	}
	double MassTerm(BasisFunction* phi1, BasisFunction* phi2)
	{
		return this->_massTerms(phi1->LocalNumber, phi2->LocalNumber);
	}

	virtual void ComputeVolumicTerm(BasisFunction* phi1, BasisFunction* phi2) = 0;
	virtual void ComputeMassTerm(BasisFunction* phi1, BasisFunction* phi2) = 0;
};