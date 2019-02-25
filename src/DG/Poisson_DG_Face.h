#pragma once
#include "../Mesh/Element.h"
#include "../FunctionalBasis/BasisFunction.h"

template <short Dim>
class Poisson_DG_Face
{
public:
	virtual double CouplingTerm(Poisson_DG_Element<Dim>* element1, BasisFunction<Dim>* phi1, Poisson_DG_Element<Dim>* element2, BasisFunction<Dim>* phi2, DiffusionPartition diffusionPartition) = 0;
	virtual double PenalizationTerm(Poisson_DG_Element<Dim>* element1, BasisFunction<Dim>* phi1, Poisson_DG_Element<Dim>* element2, BasisFunction<Dim>* phi2, double penalizationCoefficient, DiffusionPartition diffusionPartition) = 0;
};