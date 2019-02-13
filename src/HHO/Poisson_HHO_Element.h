#pragma once
#include "SourceFunction.h"
#include "BasisFunction.h"

template <short Dim>
class Poisson_HHO_Element
{
public:
	virtual double* OuterNormalVector(Face* face) = 0;
	//virtual function<double(Point)> EvalPhiOnFace(Face* face, BasisFunction<Dim>* p_phi) = 0;
	//virtual function<double*(Point)> GradPhiOnFace(Face* face, BasisFunction<Dim>* p_phi) = 0;

	virtual double InteriorTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) = 0;
	//virtual double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) = 0;
	virtual double SourceTerm(BasisFunction<Dim>* phi, SourceFunction* f) = 0;

	virtual ~Poisson_DG_Element() {}
};