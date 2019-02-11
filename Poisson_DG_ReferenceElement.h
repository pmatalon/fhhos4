#pragma once
#include <Eigen/Dense>
#include "BasisFunction.h"

template <short Dim>
class Poisson_DG_ReferenceElement
{
protected:
	Eigen::MatrixXd _volumicTerms;
	Eigen::MatrixXd _massTerms;
public:
	Poisson_DG_ReferenceElement(int nBasisFunctions) :
		_volumicTerms(nBasisFunctions, nBasisFunctions), _massTerms(nBasisFunctions, nBasisFunctions)
	{}

	double VolumicTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->_volumicTerms(phi1->LocalNumber, phi2->LocalNumber);
	}
	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->_massTerms(phi1->LocalNumber, phi2->LocalNumber);
	}

	virtual void ComputeVolumicTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) = 0;
	virtual void ComputeMassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) = 0;
};