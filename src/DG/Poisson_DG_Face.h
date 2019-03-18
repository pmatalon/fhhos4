#pragma once
#include "../Mesh/Element.h"
#include "../FunctionalBasis/BasisFunction.h"

template <short Dim>
class Poisson_DG_Face : virtual public Face<Dim>
{
public:
	Poisson_DG_Face(BigNumber number, Element<Dim>* element1, Element<Dim>* element2) : Face<Dim>(number, element1, element2) {}

	virtual double CouplingTerm(Element<Dim>* element1, BasisFunction<Dim>* phi1, Element<Dim>* element2, BasisFunction<Dim>* phi2, DiffusionPartition diffusionPartition) = 0;
	virtual double PenalizationTerm(Element<Dim>* element1, BasisFunction<Dim>* phi1, Element<Dim>* element2, BasisFunction<Dim>* phi2, double penalizationCoefficient, DiffusionPartition diffusionPartition) = 0;
};