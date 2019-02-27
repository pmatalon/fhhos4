#pragma once
#include "../Utils/SourceFunction.h"
#include "../Utils/DiffusionPartition.h"

template <short Dim>
class Poisson_DG_ReferenceElement;

template <short Dim>
class Poisson_DG_Element
{
public:
	virtual double* OuterNormalVector(Face<Dim>* face) = 0;
	virtual double DiffusionCoefficient(DiffusionPartition diffusionPartition) = 0;
	//virtual Element<Dim>* ElementOnTheOtherSideOf(Face<Dim>* face) = 0;

	virtual function<double(Point)> EvalPhiOnFace(Face<Dim>* face, BasisFunction<Dim>* phi) = 0;
	virtual function<double*(Point)> GradPhiOnFace(Face<Dim>* face, BasisFunction<Dim>* phi) = 0;

	virtual double VolumicTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2, Poisson_DG_ReferenceElement<Dim>* referenceElement, DiffusionPartition diffusionPartition) = 0;
	virtual double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2, Poisson_DG_ReferenceElement<Dim>* referenceElement) = 0;
	virtual double SourceTerm(BasisFunction<Dim>* phi, SourceFunction* f) = 0;

	virtual ~Poisson_DG_Element() {}
};