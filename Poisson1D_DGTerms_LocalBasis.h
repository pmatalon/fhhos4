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

	double VolumicTerm(Element* element, IBasisFunction1D* func1, IBasisFunction1D* func2)
	{
		Interval* interval = static_cast<Interval*>(element);
		double a = interval->A;
		double b = interval->B;

		int nQuadPoints = func1->GetDegree() + func2->GetDegree();
		if (func1->ReferenceInterval().Left == -1 && func1->ReferenceInterval().Right == 1)
		{
			// defined on [-1, 1]
			function<double(double)> functionToIntegrate = [func1, func2](double t) {
				return func1->EvalDerivative(t)*func2->EvalDerivative(t);
			};

			GaussLegendre gs(nQuadPoints);
			return 2 / (b - a) * gs.Quadrature(functionToIntegrate);
		}
		else
		{
			function<double(double)> functionToIntegrate = [func1, func2](double u) {
				//double u = 0.5 * t + 0.5;
				return func1->EvalDerivative(u)*func2->EvalDerivative(u);
			};

			return 1 / (b - a) * Utils::Integral(nQuadPoints, functionToIntegrate, 0, 1);
		}
	}

	double MassTerm(Element* element, IBasisFunction1D* func1, IBasisFunction1D* func2)
	{
		Interval* interval = static_cast<Interval*>(element);
		double a = interval->A;
		double b = interval->B;

		GaussLegendre gs(func1->GetDegree() + func2->GetDegree());

		if (func1->ReferenceInterval().Left == -1 && func1->ReferenceInterval().Right == 1)
		{
			// defined on [-1, 1]
			function<double(double)> functionToIntegrate = [func1, func2](double t) {
				return func1->Eval(t)*func2->Eval(t);
			};

			return (b - a) / 2 * gs.Quadrature(functionToIntegrate);
		}
		else
		{
			function<double(double)> functionToIntegrate = [func1, func2](double u) {
				//double u = 0.5 * t + 0.5;
				return func1->Eval(u)*func2->Eval(u);
			};

			return (b - a) * Utils::Integral(func1->GetDegree() + func2->GetDegree(), functionToIntegrate, 0, 1);
		}
	}

	double CouplingTerm(ElementInterface* interface, Element* element1, IBasisFunction1D* func1, Element* element2, IBasisFunction1D* func2)
	{
		/*if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;*/
		Interval* interval1 = static_cast<Interval*>(element1);
		Interval* interval2 = static_cast<Interval*>(element2);

		return MeanDerivative(interval1, func1, interface) * Jump(interval2, func2, interface) + MeanDerivative(interval2, func2, interface) * Jump(interval1, func1, interface);
	}

	double PenalizationTerm(ElementInterface* point, Element* element1, IBasisFunction1D* func1, Element* element2, IBasisFunction1D* func2, double penalizationCoefficient)
	{
		/*if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;*/
		Interval* interval1 = static_cast<Interval*>(element1);
		Interval* interval2 = static_cast<Interval*>(element2);

		return penalizationCoefficient * Jump(interval1, func1, point) * Jump(interval2, func2, point);
	}

	double RightHandSide(Element* element, IBasisFunction1D* func)
	{
		Interval* interval = static_cast<Interval*>(element);
		double a = interval->A;
		double b = interval->B;

		if (func->ReferenceInterval().Left == -1 && func->ReferenceInterval().Right == 1)
		{
			GaussLegendre gs;

			function<double(double)> sourceTimesBasisFunction = [this, func, a, b](double t) {
				return this->_sourceFunction((b - a) / 2 * t + (a + b) / 2) * func->Eval(t);
			};

			return (b - a) / 2 * gs.Quadrature(sourceTimesBasisFunction);
		}
		else
		{
			function<double(double)> sourceTimesBasisFunction = [this, func, a, b](double u) {
				return this->_sourceFunction((b - a) * u + a) * func->Eval(u);
			};

			return  (b - a) * Utils::Integral(sourceTimesBasisFunction, 0, 1);
		}
	}

	double MeanDerivative(Interval* element, IBasisFunction1D* func, ElementInterface* interface)
	{
		double t = interface == element->Left ? func->ReferenceInterval().Left : func->ReferenceInterval().Right; // t in [-1, 1]
		double a = element->A;
		double b = element->B;

		double jacobian = 0;
		if (func->ReferenceInterval().Left == -1 && func->ReferenceInterval().Right == 1)
			jacobian = 2 / (b - a);
		else
			jacobian = 1 / (b - a);

		double meanFactor = interface->IsDomainBoundary ? 1 : 0.5;
		return meanFactor * jacobian * func->EvalDerivative(t);
	}

	double Jump(Interval* element, IBasisFunction1D* func, ElementInterface* interface)
	{
		double t = interface == element->Left ? func->ReferenceInterval().Left : func->ReferenceInterval().Right; // t in [-1, 1]
		int factor = interface == element->Left ? 1 : -1;
		return factor * (func->Eval(t));
	}
};