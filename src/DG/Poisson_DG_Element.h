#pragma once
#include "../Utils/SourceFunction.h"
#include "../Utils/DiffusionPartition.h"

template <int Dim>
class Poisson_DG_ReferenceElement;

template <int Dim>
class Poisson_DG_Element : virtual public Element<Dim>
{
public:
	Poisson_DG_Element(BigNumber number) : Element<Dim>(number) {}

	virtual double VolumicTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2, Poisson_DG_ReferenceElement<Dim>* referenceElement, DiffusionPartition diffusionPartition) = 0;
	virtual double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2, Poisson_DG_ReferenceElement<Dim>* referenceElement) = 0;
	virtual double SourceTerm(BasisFunction<Dim>* phi, SourceFunction* f) = 0;

	virtual ~Poisson_DG_Element() {}
};