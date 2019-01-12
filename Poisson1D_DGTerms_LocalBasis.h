#pragma once
#include <functional>
#include <math.h>
#include "IBasisFunction.h"
#include "IPoisson_DGTerms.h"
#include "Utils.h"
#include "Element.h"
using namespace std;

class Poisson1D_DGTerms_LocalBasis : public IPoisson_DGTerms<IBasisFunction1D>
{
protected:
	function<double(double)> _sourceFunction;

public:
	Poisson1D_DGTerms_LocalBasis(function<double(double)> sourceFunction)
	{
		this->_sourceFunction = sourceFunction;
	}

	bool IsGlobalBasis() { return false; }

	double VolumicTerm(Element* element, IBasisFunction1D* phi1, IBasisFunction1D* phi2)
	{
		Interval* interval = static_cast<Interval*>(element);
		double h = interval->B - interval->A;

		RefInterval refInterval = phi1->ReferenceInterval();

		function<double(double)> functionToIntegrate = [phi1, phi2](double t) {
			return InnerProduct(phi1->Grad(t), phi2->Grad(t));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		double factor = refInterval.Length / h;
		return factor * Utils::Integral(nQuadPoints, functionToIntegrate, refInterval);
	}

	double MassTerm(Element* element, IBasisFunction1D* phi1, IBasisFunction1D* phi2)
	{
		Interval* interval = static_cast<Interval*>(element);
		double h = interval->B - interval->A;

		RefInterval refInterval = phi1->ReferenceInterval();

		function<double(double)> functionToIntegrate = [phi1, phi2](double t) {
			return phi1->Eval(t)*phi2->Eval(t);
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		double factor = h / refInterval.Length;
		return factor * Utils::Integral(nQuadPoints, functionToIntegrate, refInterval);
	}

	double CouplingTerm(ElementInterface* interface, Element* element1, IBasisFunction1D* phi1, Element* element2, IBasisFunction1D* phi2)
	{
		assert(interface->IsBetween(element1, element2));
		Interval* interval1 = static_cast<Interval*>(element1);
		Interval* interval2 = static_cast<Interval*>(element2);

		return MeanDerivative(interval1, phi1, interface) * Jump(interval2, phi2, interface) + MeanDerivative(interval2, phi2, interface) * Jump(interval1, phi1, interface);
	}

	double PenalizationTerm(ElementInterface* point, Element* element1, IBasisFunction1D* phi1, Element* element2, IBasisFunction1D* phi2, double penalizationCoefficient)
	{
		assert(point->IsBetween(element1, element2));
		Interval* interval1 = static_cast<Interval*>(element1);
		Interval* interval2 = static_cast<Interval*>(element2);

		return penalizationCoefficient * Jump(interval1, phi1, point) * Jump(interval2, phi2, point);
	}

	double RightHandSide(Element* element, IBasisFunction1D* phi)
	{
		Interval* interval = static_cast<Interval*>(element);
		double a = interval->A;
		double b = interval->B;

		RefInterval refInterval = phi->ReferenceInterval();

		function<double(double)> sourceTimesBasisFunction = NULL;
		if (refInterval.Left == -1 && refInterval.Right == 1)
		{
			sourceTimesBasisFunction = [this, phi, a, b](double t) {
				return this->_sourceFunction((b - a) / 2 * t + (a + b) / 2) * phi->Eval(t);
			};
		}
		else
		{
			sourceTimesBasisFunction = [this, phi, a, b](double u) {
				return this->_sourceFunction((b - a) * u + a) * phi->Eval(u);
			};
		}

		double factor = (b - a) / refInterval.Length;
		return  factor * Utils::Integral(sourceTimesBasisFunction, refInterval);
	}

	double MeanDerivative(Interval* element, IBasisFunction1D* phi, ElementInterface* interface)
	{
		RefInterval refInterval = phi->ReferenceInterval();
		double t = interface == element->Left ? refInterval.Left : refInterval.Right; // t in [-1, 1]
		double h = element->B - element->A;

		double meanFactor = interface->IsDomainBoundary ? 1 : 0.5;
		double jacobian = refInterval.Length / h;
		return meanFactor * jacobian * phi->EvalDerivative(t);
	}

	double Jump(Interval* element, IBasisFunction1D* phi, ElementInterface* interface)
	{
		double t = interface == element->Left ? phi->ReferenceInterval().Left : phi->ReferenceInterval().Right; // t in [-1, 1]
		int factor = interface == element->Left ? 1 : -1;
		return factor * (phi->Eval(t));
	}

private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0];
	}
};