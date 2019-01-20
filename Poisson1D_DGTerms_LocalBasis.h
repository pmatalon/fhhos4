#pragma once
#include <functional>
#include <math.h>
#include "IBasisFunction.h"
#include "IPoisson_DGTerms.h"
#include "Utils.h"
#include "Element.h"
#include "Poisson_DG_ReferenceInterval.h"
using namespace std;

class Poisson1D_DGTerms_LocalBasis : public IPoisson_DGTerms
{
protected:
	//function<double(double)> _sourceFunction;

public:
	Poisson1D_DGTerms_LocalBasis(function<double(double)> sourceFunction, FunctionalBasis1D* basis) 
		: IPoisson_DGTerms(new SourceFunction1D(sourceFunction))
	{
		//this->_sourceFunction = sourceFunction;
		Poisson_DG_ReferenceElement* refInterval = new Poisson_DG_ReferenceInterval(basis->NumberOfLocalFunctionsInElement(NULL));
		this->ComputeVolumicTerms(basis, refInterval);
		this->ReferenceElements.insert(std::make_pair(StandardElementCode::Interval, refInterval));
	}

	bool IsGlobalBasis() { return false; }

	/*double VolumicTerm(Element* element, IBasisFunction1D* phi1, IBasisFunction1D* phi2)
	{
		Interval* interval = static_cast<Interval*>(element);
		double h = interval->B - interval->A;

		DefInterval refInterval = phi1->DefinitionInterval();

		function<double(double)> functionToIntegrate = [phi1, phi2](double t) {
			return InnerProduct(phi1->GradPhiOnFace(t), phi2->GradPhiOnFace(t));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		double factor = refInterval.Length / h;
		double result = factor * Utils::Integral(nQuadPoints, functionToIntegrate, refInterval);
		return result;
	}*/

	double MassTerm(Element* element, IBasisFunction1D* phi1, IBasisFunction1D* phi2)
	{
		Interval* interval = static_cast<Interval*>(element);
		double h = interval->B - interval->A;

		function<double(double)> functionToIntegrate = [phi1, phi2](double t) {
			return phi1->Eval(t)*phi2->Eval(t);
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		double factor = h / 2;
		return factor * Utils::Integral(nQuadPoints, functionToIntegrate, -1,1);
	}

	/*double CouplingTerm(Face* face, Element* element1, IBasisFunction1D* phi1, Element* element2, IBasisFunction1D* phi2)
	{
		assert(face->IsBetween(element1, element2));
		Interval* interval1 = static_cast<Interval*>(element1);
		Interval* interval2 = static_cast<Interval*>(element2);

		return MeanDerivative(interval1, phi1, face) * Jump(interval2, phi2, face) + MeanDerivative(interval2, phi2, face) * Jump(interval1, phi1, face);
	}*/

	/*double PenalizationTerm(Face* point, Element* element1, IBasisFunction1D* phi1, Element* element2, IBasisFunction1D* phi2, double penalizationCoefficient) override
	{
		assert(point->IsBetween(element1, element2));
		Interval* interval1 = static_cast<Interval*>(element1);
		Interval* interval2 = static_cast<Interval*>(element2);

		return penalizationCoefficient * Jump(interval1, phi1, point) * Jump(interval2, phi2, point);
	}*/

	/*double RightHandSide(Element* element, IBasisFunction1D* phi)
	{
		Interval* interval = static_cast<Interval*>(element);
		double a = interval->A;
		double b = interval->B;

		function<double(double)> sourceTimesBasisFunction = [this, phi, a, b](double t) {
			return this->_sourceFunction((b - a) / 2 * t + (a + b) / 2) * phi->Eval(t);
		};

		double factor = (b - a) / 2;
		return  factor * Utils::Integral(sourceTimesBasisFunction, -1,1);
	}*/

	/*double MeanDerivative(Interval* element, IBasisFunction1D* phi, Face* interface)
	{
		double t = interface == element->Left ? -1 : 1; // t in [-1, 1]
		double h = element->B - element->A;

		double meanFactor = interface->IsDomainBoundary ? 1 : 0.5;
		double jacobian = 2 / h;
		return meanFactor * jacobian * phi->EvalDerivative(t);
	}

	double Jump(Interval* element, IBasisFunction1D* phi, Face* point)
	{
		double t = point == element->Left ? -1 : 1; // t in [-1, 1]
		int factor = point == element->Left ? 1 : -1;
		return factor * (phi->Eval(t));
	}*/

private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0];
	}
};