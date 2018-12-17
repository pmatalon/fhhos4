#pragma once
#include <map>
#include <functional>
#include <math.h>
#include "FunctionalBasisWithObjects.h"
#include "CartesianGrid2D.h"
#include "IBasisFunction2D.h"
#include "Utils.h"
#include "Element.h"
#include "ElementInterface.h"
#include <cstdio>
#include <assert.h>
using namespace std;

class FunctionalBasis2D : public FunctionalBasisWithObjects<IBasisFunction2D>
{
protected:
	int _penalizationCoefficient;
	function<double(double, double)> _sourceFunction;

	map<int, IBasisFunction2D*> _localFunctions;

public:
	FunctionalBasis2D(int penalizationCoefficient, function<double(double, double)> sourceFunction)
	{
		this->_penalizationCoefficient = penalizationCoefficient;
		this->_sourceFunction = sourceFunction;
	}

	int NumberOfLocalFunctionsInElement(Element* element)
	{
		return static_cast<int>(this->_localFunctions.size());
	}

	IBasisFunction2D* GetLocalBasisFunction(Element* element, int localFunctionNumber)
	{
		return this->_localFunctions[localFunctionNumber];
	}

	BigNumber GlobalFunctionNumber(Element* element, int localFunctionNumber)
	{
		return element->Number * static_cast<int>(this->_localFunctions.size()) + localFunctionNumber + 1; // +1 so that the numbers start at 1
	}

	double VolumicTerm(Element* element, IBasisFunction2D* phi1, IBasisFunction2D* phi2)
	{
		Square* square = static_cast<Square*>(element);
		return this->VolumicTerm(square, phi1, phi2);
	}

	double VolumicTerm(Square* element, IBasisFunction2D* phi1, IBasisFunction2D* phi2)
	{
		double h = element->Width;

		function<double(double, double)> functionToIntegrate = [phi1, phi2, h](double t, double u) {
			return phi1->EvalGradX(t, u)*phi2->EvalGradX(t, u) + phi1->EvalGradY(t, u)*phi2->EvalGradY(t, u);
		};

		GaussLegendre gs(phi1->GetDegree() + phi2->GetDegree());
		//return 4 / pow(h, 2) * gs.Quadrature(functionToIntegrate);
		//return pow(h, 2) / 4 * gs.Quadrature(functionToIntegrate);
		return gs.Quadrature(functionToIntegrate);
	}

	double CouplingTerm(ElementInterface* interface, Element* element1, IBasisFunction2D* phi1, Element* element2, IBasisFunction2D* phi2)
	{
		Square* square1 = static_cast<Square*>(element1);
		Square* square2 = static_cast<Square*>(element2);
		return this->CouplingTerm(interface, square1, phi1, square2, phi2);
	}

	double CouplingTerm(ElementInterface* interface, Square* element1, IBasisFunction2D* phi1, Square* element2, IBasisFunction2D* phi2)
	{
		if (!interface->IsBetween(element1, element2))
			return 0;

		double h = element1->Width;

		auto n1 = element1->OuterNormalVector(interface);
		auto n2 = element2->OuterNormalVector(interface);
		//auto n2 = n1;

		assert(n1[0] != 0 || n1[1] != 0);
		assert(n2[0] != 0 || n2[1] != 0);
		assert((n1[0] == n2[0] && n1[1] == n2[1]) || (n1[0] == -n2[0] && n1[1] == -n2[1]));

		double meanFactor = interface->IsDomainBoundary ? 1 : 0.5;

		Element2DInterface* interf = (Element2DInterface*)interface;

		/*function<double(double, double)> functionToIntegrate = [meanFactor, n1, n2, phi1, phi2, h, element1, element2](double t, double u) {
			double meanGradPhi1_scal_jumpPhi2 = meanFactor * (phi1->EvalGradX(t, u) * n2[0] + phi1->EvalGradY(t, u) * n2[1]) * phi2->Eval(t, u);
			double meanGradPhi2_scal_jumpPhi1 = meanFactor * (phi2->EvalGradX(t, u) * n1[0] + phi2->EvalGradY(t, u) * n1[1]) * phi1->Eval(t, u);
			return 2/h*(meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1);
		};*/


		GaussLegendre gs(phi1->GetDegree() + phi2->GetDegree() + 1);
		std::function<double(double)> functionToIntegrate;
		if (interf->IsVertical())
		{
			double t1 = interface == element1->EastInterface ? 1 : -1;
			double t2 = interface == element2->EastInterface ? 1 : -1;
			functionToIntegrate = [meanFactor, n1, n2, phi1, phi2, h, t1, t2](double u) {
				double meanGradPhi1_scal_jumpPhi2 = meanFactor * (phi1->EvalGradX(t1, u) * n2[0] + phi1->EvalGradY(t1, u) * n2[1]) * phi2->Eval(t2, u);
				double meanGradPhi2_scal_jumpPhi1 = meanFactor * (phi2->EvalGradX(t2, u) * n1[0] + phi2->EvalGradY(t2, u) * n1[1]) * phi1->Eval(t1, u);
				return 2 / h * (meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1);
			};
		}
		else if (interf->IsHorizontal())
		{
			double u1 = interface == element1->NorthInterface ? 1 : -1;
			double u2 = interface == element2->NorthInterface ? 1 : -1;
			functionToIntegrate = [meanFactor, n1, n2, phi1, phi2, h, u1, u2](double t) {
				double meanGradPhi1_scal_jumpPhi2 = meanFactor * (phi1->EvalGradX(t, u1) * n2[0] + phi1->EvalGradY(t, u1) * n2[1]) * phi2->Eval(t, u2);
				double meanGradPhi2_scal_jumpPhi1 = meanFactor * (phi2->EvalGradX(t, u2) * n1[0] + phi2->EvalGradY(t, u2) * n1[1]) * phi1->Eval(t, u1);
				return 2 / h * (meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1);
			};
		}
		else
			return 0;

		return -h / 2 * gs.Quadrature(functionToIntegrate);
		//return -pow(h, 2) / 4 * gs.Quadrature(func1D);
		//return -2/h* gs.Quadrature(func1D);
		//return -gs.Quadrature(func1D);
	}

	double PenalizationTerm(ElementInterface* interface, Element* element1, IBasisFunction2D* phi1, Element* element2, IBasisFunction2D* phi2)
	{
		Square* square1 = static_cast<Square*>(element1);
		Square* square2 = static_cast<Square*>(element2);
		return this->PenalizationTerm(interface, square1, phi1, square2, phi2);
	}

	double PenalizationTerm(ElementInterface* interface, Square* element1, IBasisFunction2D* phi1, Square* element2, IBasisFunction2D* phi2)
	{
		//if (!interface->IsBetween(element1, element2))
		//	return 0;
		double h = element1->Width;
		auto n1 = element1->OuterNormalVector(interface);
		auto n2 = element2->OuterNormalVector(interface);

		assert(n1[0] != 0 || n1[1] != 0);
		assert(n2[0] != 0 || n2[1] != 0);
		assert((n1[0] == n2[0] && n1[1] == n2[1]) || (n1[0] == -n2[0] && n1[1] == -n2[1]));

		Element2DInterface* interf = (Element2DInterface*)interface;

		/*function<double(double, double)> functionToIntegrate = [phi1, phi2, n1, n2, interf, element1, element2](double t, double u) {
			return (n1[0] * n2[0] + n1[1] * n2[1]) * phi1->Eval(t, u) * phi2->Eval(t, u);
		};*/

		double integralJump1ScalarJump2 = 0;

		
		GaussLegendre gs(phi1->GetDegree() + phi2->GetDegree() + 2);
		//std::function<double(double)> func1D;
		function<double(double)> functionToIntegrate;

		assert(interf->IsVertical() || interf->IsHorizontal());

		if (interf->IsVertical())
		{
			double t1 = interface == element1->EastInterface ? 1 : -1;
			double t2 = interface == element2->EastInterface ? 1 : -1;
			/*func1D = [functionToIntegrate, t](double u) {
				return functionToIntegrate(t, u);
			};*/

			functionToIntegrate = [phi1, phi2, n1, n2, t1, t2](double u) {
				return (n1[0] * n2[0] + n1[1] * n2[1]) * phi1->Eval(t1, u) * phi2->Eval(t2, u);
			};
		}
		else if (interf->IsHorizontal())
		{
			double u1 = interface == element1->NorthInterface ? 1 : -1;
			double u2 = interface == element2->NorthInterface ? 1 : -1;
			functionToIntegrate = [phi1, phi2, n1, n2, u1, u2](double t) {
				return (n1[0] * n2[0] + n1[1] * n2[1]) * phi1->Eval(t, u1) * phi2->Eval(t, u2);
			};
		}
		else
			return 0;

		integralJump1ScalarJump2 = h / 2 * gs.Quadrature(functionToIntegrate);
		//integralJump1ScalarJump2 = pow(h, 2) / 4 * gs.Quadrature(func1D);
		//integralJump1ScalarJump2 = 2 / h * gs.Quadrature(func1D);
		//integralJump1ScalarJump2 = gs.Quadrature(func1D);
		return this->_penalizationCoefficient * integralJump1ScalarJump2;
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

		function<double(double, double)> sourceTimesBasisFunction = [this, phi, x1, x2, y1, y2](double t, double u) {
			assert((x2 - x1) / 2 * t + (x2 + x1) / 2 >= x1 && (x2 - x1) / 2 * t + (x2 + x1) / 2 <= x2);
			assert((y2 - y1) / 2 * u + (y2 + y1) / 2 >= y1 && (y2 - y1) / 2 * u + (y2 + y1) / 2 <= y2);
			return this->_sourceFunction((x2 - x1) / 2 * t + (x2 + x1) / 2, (y2 - y1) / 2 * u + (y2 + y1) / 2) * phi->Eval(t, u);
		};

		cout << "abs = " << abs((x2 - x1) * (y2 - y1) / 4 - pow(element->Width, 2) / 4) << endl;
		assert(abs((x2 - x1) * (y2 - y1) / 4 - pow(element->Width, 2)/4) < pow(10, -10));

		GaussLegendre gs;
		return (x2 - x1) * (y2 - y1) / 4 * gs.Quadrature(sourceTimesBasisFunction);
	}
};