#pragma once
#include "Face.h"
#include "Poisson_DG_Face.h"
#include "Interval.h"

class Face1D : public Face, public Poisson_DG_Face
{
public:
	double X;

	/*Face1D(BigNumber number, double x, Element* element1, Element* element2) : Face(number, element1, element2)
	{	
		this->X = x;
	}

	Face1D(BigNumber number, double x, Element* element1) : Face(number, element1)
	{	
		this->X = x;
	}*/

	Face1D(BigNumber number, double x) : Face(number, NULL, NULL)
	{
		this->X = x;
		this->IsDomainBoundary = false;
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double CouplingTerm(Poisson_DG_Element* element1, BasisFunction* p_phi1, Poisson_DG_Element* element2, BasisFunction* p_phi2)
	{
		IBasisFunction1D* phi1 = static_cast<IBasisFunction1D*>(p_phi1);
		IBasisFunction1D* phi2 = static_cast<IBasisFunction1D*>(p_phi2);

		Interval* interval1 = dynamic_cast<Interval*>(element1);
		Interval* interval2 = dynamic_cast<Interval*>(element2);

		return MeanDerivative(interval1, phi1) * Jump(interval2, phi2) + MeanDerivative(interval2, phi2) * Jump(interval1, phi1);
	}

	double PenalizationTerm(Poisson_DG_Element* element1, BasisFunction* p_phi1, Poisson_DG_Element* element2, BasisFunction* p_phi2, double penalizationCoefficient)
	{
		IBasisFunction1D* phi1 = static_cast<IBasisFunction1D*>(p_phi1);
		IBasisFunction1D* phi2 = static_cast<IBasisFunction1D*>(p_phi2);

		Interval* interval1 = static_cast<Interval*>(element1);
		Interval* interval2 = static_cast<Interval*>(element2);

		return penalizationCoefficient * Jump(interval1, phi1) * Jump(interval2, phi2);
	}

	double MeanDerivative(Interval* element, IBasisFunction1D* phi)
	{
		double t = this == element->Left ? -1 : 1; // t in [-1, 1]
		double h = element->B - element->A;

		double meanFactor = this->IsDomainBoundary ? 1 : 0.5;
		double jacobian = 2 / h;
		return meanFactor * jacobian * phi->EvalDerivative(t);
	}

	double Jump(Interval* element, IBasisFunction1D* phi)
	{
		double t = this == element->Left ? -1 : 1; // t in [-1, 1]
		int factor = this == element->Left ? 1 : -1;
		return factor * (phi->Eval(t));
	}
};