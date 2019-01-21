#pragma once
#include "SourceFunction.h"

class Poisson_DG_ReferenceElement;

class Poisson_DG_Element
{
public:
	virtual double* OuterNormalVector(Face* face) = 0;
	virtual function<double(Point)> EvalPhiOnFace(Face* face, BasisFunction* p_phi) = 0;
	virtual function<double*(Point)> GradPhiOnFace(Face* face, BasisFunction* p_phi) = 0;

	virtual double VolumicTerm(BasisFunction* phi1, BasisFunction* phi2, Poisson_DG_ReferenceElement* referenceElement) = 0;
	virtual double MassTerm(BasisFunction* phi1, BasisFunction* phi2, Poisson_DG_ReferenceElement* referenceElement) = 0;
	virtual double SourceTerm(BasisFunction* phi, SourceFunction* f) = 0;

	virtual ~Poisson_DG_Element() {}
};