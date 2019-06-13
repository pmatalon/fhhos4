#pragma once
#include "../Face.h"
#include "../../DG/Poisson_DG_Face.h"
#include "Interval.h"

class InterfacePoint : virtual public Face<1>, public Poisson_DG_Face<1>
{
public:
	double X;

	InterfacePoint(BigNumber number, double x) : Face(number, NULL, NULL), Poisson_DG_Face(number, NULL, NULL)
	{
		this->X = x;
		this->IsDomainBoundary = false;
	}

	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	double GetDiameter()
	{
		if (this->Element1 != NULL)
			return this->Element1->GetDiameter();
		else
			return this->Element2->GetDiameter();
	}

	double Measure()
	{
		return 0;
	}

	virtual double MassTerm(BasisFunction<0>* facePhi, Element<1>* element, BasisFunction<1>* reconstructPhi)
	{
		return 0;
	}

	RefPoint ConvertToReference(DomPoint domainPoint)
	{
		return RefPoint(0);
	}

	DomPoint ConvertToDomain(RefPoint referenceElementPoint)
	{
		return DomPoint(0);
	}

	vector<RefPoint> GetNodalPoints(FunctionalBasis<0>* basis)
	{
		return vector<RefPoint> {RefPoint(0)};
	}

	double ComputeIntegral(function<double(RefPoint)> func, int polynomialDegree)
	{
		assert(false);
	}

	double ComputeIntegral(function<double(RefPoint)> func)
	{
		assert(false);
	}

	//---------------------------------------------------------------//
	//                 Poisson_DG_Face implementation                //
	//---------------------------------------------------------------//

	double CouplingTerm(Element<1>* element1, BasisFunction<1>* p_phi1, Element<1>* element2, BasisFunction<1>* p_phi2, DiffusionPartition diffusionPartition) override
	{
		IBasisFunction1D* phi1 = static_cast<IBasisFunction1D*>(p_phi1);
		IBasisFunction1D* phi2 = static_cast<IBasisFunction1D*>(p_phi2);

		Interval* interval1 = dynamic_cast<Interval*>(element1);
		Interval* interval2 = dynamic_cast<Interval*>(element2);

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

	double PenalizationTerm(Element<1>* element1, BasisFunction<1>* p_phi1, Element<1>* element2, BasisFunction<1>* p_phi2, double penalizationCoefficient, DiffusionPartition diffusionPartition) override
	{
		IBasisFunction1D* phi1 = static_cast<IBasisFunction1D*>(p_phi1);
		IBasisFunction1D* phi2 = static_cast<IBasisFunction1D*>(p_phi2);

		Interval* interval1 = dynamic_cast<Interval*>(element1);
		Interval* interval2 = dynamic_cast<Interval*>(element2);

		double diffusionDependantCoefficient = element1->DiffusionCoefficient(diffusionPartition);
		if (!this->IsDomainBoundary)
		{
			double k1 = this->Element1->DiffusionCoefficient(diffusionPartition);
			double k2 = this->Element2->DiffusionCoefficient(diffusionPartition);
			diffusionDependantCoefficient = 2 * k1*k2 / (k1 + k2);
		}

		return diffusionDependantCoefficient * penalizationCoefficient * Jump(interval1, phi1) * Jump(interval2, phi2);
	}

	double MeanDerivative(Interval* element, IBasisFunction1D* phi)
	{
		double t = this == element->Left ? -1 : 1; // t in [-1, 1]
		double h = element->WidthX;
		return 2 / h * phi->EvalDerivative(t);
	}

	double Jump(Interval* element, IBasisFunction1D* phi)
	{
		double t = this == element->Left ? -1 : 1; // t in [-1, 1]
		int normalFactor = this == element->Left ? 1 : -1;
		return normalFactor * (phi->Eval(t));
	}
};