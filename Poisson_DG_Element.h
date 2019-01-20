#pragma once
#include "Point.h"

class Poisson_DG_ReferenceElement;

class Poisson_DG_Element
{
public:
	virtual double VolumicTerm(BasisFunction* phi1, BasisFunction* phi2, Poisson_DG_ReferenceElement* referenceElement) = 0;

	virtual double* OuterNormalVector(Face* face) = 0;
	virtual function<double(Point)> EvalPhiOnFace(Face* face, BasisFunction* p_phi) = 0;
	virtual function<double*(Point)> GradPhiOnFace(Face* face, BasisFunction* p_phi) = 0;

	//virtual double RightHandSide(BasisFunction* phi1, Poisson_DG_ReferenceElement* referenceElement) = 0;
	virtual ~Poisson_DG_Element() {}
};