#pragma once
#include "../Utils/SourceFunction.h"
#include "../FunctionalBasis/BasisFunction.h"

template <short Dim>
class Poisson_HHO_Element
{
public:
	virtual double* OuterNormalVector(Face<Dim>* face) = 0;
	//virtual function<double(Point)> EvalPhiOnFace(Face* face, BasisFunction<Dim>* p_phi) = 0;
	//virtual function<double*(Point)> GradPhiOnFace(Face* face, BasisFunction<Dim>* p_phi) = 0;

	//virtual double InteriorTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) = 0;
	virtual double S(BasisFunction<Dim>* phiReconstruct1, BasisFunction<Dim>* phiReconstruct2) = 0;
	virtual double Bt(BasisFunction<Dim>* phiReconstruct, BasisFunction<Dim>* phiElement) = 0;
	virtual double Bf(BasisFunction<Dim>* phiReconstruct, BasisFunction<Dim-1>* phiFace) = 0;
	//virtual double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) = 0;
	virtual double SourceTerm(BasisFunction<Dim>* phi, SourceFunction* f) = 0;

	virtual ~Poisson_HHO_Element() {}
};