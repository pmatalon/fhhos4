#pragma once
#include "Face.h"
#include "../DG/Poisson_DG_Face.h"
#include "Interval.h"

class PointFace : public Face<1>, public Poisson_DG_Face<1>
{
public:
	double X;

	PointFace(BigNumber number, double x) : Face(number, NULL, NULL)
	{
		this->X = x;
		this->IsDomainBoundary = false;
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double CouplingTerm(Element<1>* element1, BasisFunction<1>* p_phi1, Element<1>* element2, BasisFunction<1>* p_phi2, DiffusionPartition diffusionPartition)
	{
		IBasisFunction1D* phi1 = static_cast<IBasisFunction1D*>(p_phi1);
		IBasisFunction1D* phi2 = static_cast<IBasisFunction1D*>(p_phi2);

		Interval* interval1 = static_cast<Interval*>(element1);
		Interval* interval2 = static_cast<Interval*>(element2);

		double k1 = element1->DiffusionCoefficient(diffusionPartition);
		double k2 = element2->DiffusionCoefficient(diffusionPartition);

		double weight1 = 1;
		double weight2 = 1;
		if (!this->IsDomainBoundary)
		{
			Element<1>* elementOnTheOtherSide1 = interval1->ElementOnTheOtherSideOf(this);
			Element<1>* elementOnTheOtherSide2 = interval2->ElementOnTheOtherSideOf(this);
			double l1 = k1;
			double l2 = elementOnTheOtherSide1->DiffusionCoefficient(diffusionPartition);
			weight1 = elementOnTheOtherSide1->DiffusionCoefficient(diffusionPartition) / (l1 + l2);
			weight2 = elementOnTheOtherSide2->DiffusionCoefficient(diffusionPartition) / (l1 + l2);
		}

		return weight1 * k1 * MeanDerivative(interval1, phi1) * Jump(interval2, phi2) + weight2 * k2 * MeanDerivative(interval2, phi2) * Jump(interval1, phi1);
	}

	double PenalizationTerm(Poisson_DG_Element<1>* element1, BasisFunction<1>* p_phi1, Poisson_DG_Element<1>* element2, BasisFunction<1>* p_phi2, double penalizationCoefficient, DiffusionPartition diffusionPartition)
	{
		IBasisFunction1D* phi1 = static_cast<IBasisFunction1D*>(p_phi1);
		IBasisFunction1D* phi2 = static_cast<IBasisFunction1D*>(p_phi2);

		Interval* interval1 = static_cast<Interval*>(element1);
		Interval* interval2 = static_cast<Interval*>(element2);

		double diffusionDependantCoefficient = element1->DiffusionCoefficient(diffusionPartition);
		if (!this->IsDomainBoundary)
		{
			double k1 = this->Element1->DiffusionCoefficient(diffusionPartition);
			double k2 = this->Element2->DiffusionCoefficient(diffusionPartition);
			diffusionDependantCoefficient = 2 * k1*k2 / (k1 + k2);
		}
		//double h = interval1->B - interval1->A;
		return diffusionDependantCoefficient * penalizationCoefficient * Jump(interval1, phi1) * Jump(interval2, phi2);
	}

	double MeanDerivative(Interval* element, IBasisFunction1D* phi)
	{
		double t = this == element->Left ? -1 : 1; // t in [-1, 1]
		double h = element->B - element->A;
		return 2 / h * phi->EvalDerivative(t);
	}

	double Jump(Interval* element, IBasisFunction1D* phi)
	{
		double t = this == element->Left ? -1 : 1; // t in [-1, 1]
		int normalFactor = this == element->Left ? 1 : -1;
		return normalFactor * (phi->Eval(t));
	}
};