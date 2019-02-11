#pragma once
#include "SourceFunction.h"

template <short Dim>
class Poisson_DG_ReferenceElement;

template <short Dim>
class Poisson_DG_Element
{
public:
	virtual double* OuterNormalVector(Face* face) = 0;
	virtual function<double(Point)> EvalPhiOnFace(Face* face, BasisFunction<Dim>* p_phi) = 0;
	virtual function<double*(Point)> GradPhiOnFace(Face* face, BasisFunction<Dim>* p_phi) = 0;

	virtual double VolumicTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2, Poisson_DG_ReferenceElement<Dim>* referenceElement) = 0;
	virtual double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2, Poisson_DG_ReferenceElement<Dim>* referenceElement) = 0;
	virtual double SourceTerm(BasisFunction<Dim>* phi, SourceFunction* f) = 0;

	virtual ~Poisson_DG_Element() {}
};