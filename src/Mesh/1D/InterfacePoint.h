#pragma once
#include "Interval.h"
#include "../CartesianFace.h"

class InterfacePoint : public CartesianFace<1>
{
public:
	Vertex* V;

	InterfacePoint(BigNumber number, Vertex* v) : 
		Face(number, NULL, NULL), 
		CartesianFace<1>(number, v, 0, NULL, NULL, CartesianShapeOrientation::None)
	{
		this->V = v;
		this->IsDomainBoundary = false;
	}

	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	double Diameter() override
	{
		if (this->Element1 != NULL)
			return this->Element1->Diameter();
		else
			return this->Element2->Diameter();
	}

	double Measure() override
	{
		return 0;
	}

	Face<1>* CreateSameGeometricFace(BigNumber number, Element<1>* element1)
	{
		return new InterfacePoint(number, this->V);
	}

	void ExportFaceToMatlab(FILE* file)
	{
		fprintf(file, "%llu %.17g 0 0 0 %d\n", this->Number, this->V->X, this->IsDomainBoundary);
	}

	//---------------------------------------------------------------//
	//                 Poisson_DG_Face implementation                //
	//---------------------------------------------------------------//

	double CouplingTerm(Element<1>* element1, BasisFunction<1>* p_phi1, Element<1>* element2, BasisFunction<1>* p_phi2) override
	{
		IBasisFunction1D* phi1 = static_cast<IBasisFunction1D*>(p_phi1);
		IBasisFunction1D* phi2 = static_cast<IBasisFunction1D*>(p_phi2);

		Interval* interval1 = dynamic_cast<Interval*>(element1);
		Interval* interval2 = dynamic_cast<Interval*>(element2);

		double k1 = element1->Kappa;
		double k2 = element2->Kappa;

		double weight1 = 1;
		double weight2 = 1;
		if (!this->IsDomainBoundary)
		{
			Element<1>* elementOnTheOtherSide1 = interval1->ElementOnTheOtherSideOf(this);
			Element<1>* elementOnTheOtherSide2 = interval2->ElementOnTheOtherSideOf(this);
			double l1 = k1;
			double l2 = elementOnTheOtherSide1->Kappa;
			weight1 = elementOnTheOtherSide1->Kappa / (l1 + l2);
			weight2 = elementOnTheOtherSide2->Kappa / (l1 + l2);
		}

		return weight1 * k1 * MeanDerivative(interval1, phi1) * Jump(interval2, phi2) + weight2 * k2 * MeanDerivative(interval2, phi2) * Jump(interval1, phi1);
	}

	double PenalizationTerm(Element<1>* element1, BasisFunction<1>* p_phi1, Element<1>* element2, BasisFunction<1>* p_phi2, double penalizationCoefficient) override
	{
		IBasisFunction1D* phi1 = static_cast<IBasisFunction1D*>(p_phi1);
		IBasisFunction1D* phi2 = static_cast<IBasisFunction1D*>(p_phi2);

		Interval* interval1 = dynamic_cast<Interval*>(element1);
		Interval* interval2 = dynamic_cast<Interval*>(element2);

		double diffusionDependantCoefficient = element1->Kappa;
		if (!this->IsDomainBoundary)
		{
			double k1 = this->Element1->Kappa;
			double k2 = this->Element2->Kappa;
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