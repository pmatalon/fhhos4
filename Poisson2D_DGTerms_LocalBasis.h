#pragma once
#include <functional>
#include <math.h>
#include "IPoisson_DGTerms.h"
#include "IBasisFunction.h"
#include "Utils.h"
#include "Element.h"
#include "ElementInterface.h"
#include "Square.h"
using namespace std;

class Poisson2D_DGTerms_LocalBasis : public IPoisson_DGTerms<IBasisFunction2D>
{
protected:
	function<double(double, double)> _sourceFunction;

public:
	Poisson2D_DGTerms_LocalBasis(function<double(double, double)> sourceFunction)
	{
		this->_sourceFunction = sourceFunction;
	}

	bool IsGlobalBasis() { return false; }

	double VolumicTerm(Element* element, IBasisFunction2D* phi1, IBasisFunction2D* phi2)
	{
		Square* square = static_cast<Square*>(element);
		return this->VolumicTerm(square, phi1, phi2);
	}

	double VolumicTerm(Square* element, IBasisFunction2D* phi1, IBasisFunction2D* phi2)
	{
		function<double(double, double)> functionToIntegrate = [phi1, phi2](double t, double u) {
			return InnerProduct(phi1->Grad(t, u), phi2->Grad(t, u));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		return Utils::Integral(nQuadPoints, functionToIntegrate, phi1->ReferenceInterval(), phi1->ReferenceInterval());
	}

	double CouplingTerm(ElementInterface* interface, Element* element1, IBasisFunction2D* phi1, Element* element2, IBasisFunction2D* phi2)
	{
		assert(interface->IsBetween(element1, element2));

		Square* square1 = static_cast<Square*>(element1);
		Square* square2 = static_cast<Square*>(element2);
		return this->CouplingTerm(interface, square1, phi1, square2, phi2);
	}

	double CouplingTerm(ElementInterface* interface, Square* element1, IBasisFunction2D* phi1, Square* element2, IBasisFunction2D* phi2)
	{
		RefInterval refInterval = phi1->ReferenceInterval();
		
		auto n1 = element1->OuterNormalVector(interface);
		auto n2 = element2->OuterNormalVector(interface);

		double meanFactor = interface->IsDomainBoundary ? 1 : 0.5;

		Element2DInterface* interf = (Element2DInterface*)interface;

		std::function<double(double)> functionToIntegrate;
		if (interf->IsVertical())
		{
			double t1 = interface == element1->EastInterface ? refInterval.Right : refInterval.Left;
			double t2 = interface == element2->EastInterface ? refInterval.Right : refInterval.Left;

			functionToIntegrate = [n1, n2, phi1, phi2, t1, t2](double u) {
				double meanGradPhi1_scal_jumpPhi2 = InnerProduct(phi1->Grad(t1, u), n2) * phi2->Eval(t2, u);
				double meanGradPhi2_scal_jumpPhi1 = InnerProduct(phi2->Grad(t2, u), n1) * phi1->Eval(t1, u);
				return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
			};
		}
		else if (interf->IsHorizontal())
		{
			double u1 = interface == element1->NorthInterface ? refInterval.Right : refInterval.Left;
			double u2 = interface == element2->NorthInterface ? refInterval.Right : refInterval.Left;

			functionToIntegrate = [n1, n2, phi1, phi2, u1, u2](double t) {
				double meanGradPhi1_scal_jumpPhi2 = InnerProduct(phi1->Grad(t, u1), n2) * phi2->Eval(t, u2);
				double meanGradPhi2_scal_jumpPhi1 = InnerProduct(phi2->Grad(t, u2), n1) * phi1->Eval(t, u1);
				return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
			};
		}
		else
			return 0;

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 1;
		return -meanFactor * Utils::Integral(nQuadPoints, functionToIntegrate, refInterval);
	}

	double PenalizationTerm(ElementInterface* interface, Element* element1, IBasisFunction2D* phi1, Element* element2, IBasisFunction2D* phi2, double penalizationCoefficient)
	{
		assert(interface->IsBetween(element1, element2));

		Square* square1 = static_cast<Square*>(element1);
		Square* square2 = static_cast<Square*>(element2);
		return this->PenalizationTerm(interface, square1, phi1, square2, phi2, penalizationCoefficient);
	}

	double PenalizationTerm(ElementInterface* interface, Square* element1, IBasisFunction2D* phi1, Square* element2, IBasisFunction2D* phi2, double penalizationCoefficient)
	{
		double h = element1->Width;
		auto n1 = element1->OuterNormalVector(interface);
		auto n2 = element2->OuterNormalVector(interface);

		RefInterval refInterval = phi1->ReferenceInterval();

		Element2DInterface* interf = (Element2DInterface*)interface;

		function<double(double)> functionToIntegrate;

		if (interf->IsVertical())
		{
			double t1 = interface == element1->EastInterface ? refInterval.Right : refInterval.Left;
			double t2 = interface == element2->EastInterface ? refInterval.Right : refInterval.Left;

			functionToIntegrate = [phi1, phi2, n1, n2, t1, t2](double u) {
				return InnerProduct(n1, n2) * phi1->Eval(t1, u) * phi2->Eval(t2, u);
			};
		}
		else if (interf->IsHorizontal())
		{
			double u1 = interface == element1->NorthInterface ? refInterval.Right : refInterval.Left;
			double u2 = interface == element2->NorthInterface ? refInterval.Right : refInterval.Left;

			functionToIntegrate = [phi1, phi2, n1, n2, u1, u2](double t) {
				return InnerProduct(n1, n2) * phi1->Eval(t, u1) * phi2->Eval(t, u2);
			};
		}
		else
			return 0;

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		double jacobian = h / refInterval.Length;
		double integralJump1ScalarJump2 = jacobian * Utils::Integral(nQuadPoints, functionToIntegrate, refInterval);
		return penalizationCoefficient * integralJump1ScalarJump2;
	}

	double RightHandSide(Element* element, IBasisFunction2D* phi)
	{
		Square* square = static_cast<Square*>(element);
		return this->RightHandSide(square, phi);
	}

	double RightHandSide(Square* element, IBasisFunction2D* phi)
	{
		double x1 = element->X;
		double x2 = element->X + element->Width;
		double y1 = element->Y;
		double y2 = element->Y + element->Width;

		RefInterval refInterval = phi->ReferenceInterval();

		function<double(double, double)> sourceTimesBasisFunction = NULL;
		if (refInterval.Left == -1 && refInterval.Right == 1)
		{
			sourceTimesBasisFunction = [this, phi, x1, x2, y1, y2](double t, double u) {
				return this->_sourceFunction((x2 - x1) / 2 * t + (x2 + x1) / 2, (y2 - y1) / 2 * u + (y2 + y1) / 2) * phi->Eval(t, u);
			};
		}
		else
		{
			sourceTimesBasisFunction = [this, phi, x1, x2, y1, y2](double t, double u) {
				return this->_sourceFunction((x2 - x1) * t + x1, (y2 - y1) * u + y1) * phi->Eval(t, u);
			};
		}

		double jacobian = (x2 - x1) * (y2 - y1) / pow(refInterval.Length, 2);
		return jacobian * Utils::Integral(sourceTimesBasisFunction, refInterval, refInterval);
	}

private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0] + vector1[1] * vector2[1];
	}
};