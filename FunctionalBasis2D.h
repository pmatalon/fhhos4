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

	double VolumicTerm(Element* element, IBasisFunction2D* func1, IBasisFunction2D* func2)
	{
		Square* square = static_cast<Square*>(element);
		return this->VolumicTerm(square, func1, func2);
	}

	double VolumicTerm(Square* element, IBasisFunction2D* func1, IBasisFunction2D* func2)
	{
		double h = element->Width;

		function<double(double, double)> functionToIntegrate = [func1, func2, h](double t, double u) {
			return /*4/pow(h, 2) **/ (func1->EvalGradX(t, u)*func2->EvalGradX(t, u) + func1->EvalGradY(t, u)*func2->EvalGradY(t, u));
		};

		IPolynomialFunction* poly1 = dynamic_cast<IPolynomialFunction*>(func1);
		IPolynomialFunction* poly2 = dynamic_cast<IPolynomialFunction*>(func2);

		GaussLegendre* gs = new GaussLegendre(poly1->GetDegree() + poly2->GetDegree());
		//return 4 / pow(h, 2) * gs->Quadrature(functionToIntegrate);
		//return pow(h, 2) / 4 * gs->Quadrature(functionToIntegrate);
		return gs->Quadrature(functionToIntegrate);
	}

	double CouplingTerm(ElementInterface* interface, Element* element1, IBasisFunction2D* func1, Element* element2, IBasisFunction2D* func2)
	{
		if (!interface->IsBetween(element1, element2))
			return 0;

		double h = ((Square*)element1)->Width;

		auto n1 = element1->OuterNormalVector(interface);
		auto n2 = element2->OuterNormalVector(interface);
		//auto n2 = n1;

		double meanFactor = interface->IsDomainBoundary ? 1 : 0.5;
		/*double jumpFactor = interface->IsDomainBoundary ? 1 : -1;

		// {{grad f1}}_x
		function<double(double, double)> meanGradFunc1X = [func1, meanFactor, h](double t, double u) {
			return meanFactor * 2 / h * func1->EvalGradX(t, u);
		};
		// {{grad f1}}_y
		function<double(double, double)> meanGradFunc1Y = [func1, meanFactor, h](double t, double u) {
			return meanFactor * 2 / h * func1->EvalGradY(t, u);
		};
		// {{grad f2}}_x
		function<double(double, double)> meanGradFunc2X = [func2, meanFactor, h](double t, double u) {
			return meanFactor * 2 / h * func2->EvalGradX(t, u);
		};
		// {{grad f2}}_y
		function<double(double, double)> meanGradFunc2Y = [func2, meanFactor, h](double t, double u) {
			return meanFactor * 2 / h * func2->EvalGradY(t, u);
		};
		// [[f1]]_x
		function<double(double, double)> jumpFunc1X = [func1, n1, jumpFactor](double t, double u) {
			return jumpFactor * func1->Eval(t, u)*n1[0];
		};
		// [[f1]]_y
		function<double(double, double)> jumpFunc1Y = [func1, n1, jumpFactor](double t, double u) {
			return jumpFactor * func1->Eval(t, u)*n1[1];
		};
		// [[f2]]_x
		function<double(double, double)> jumpFunc2X = [func2, n2, jumpFactor](double t, double u) {
			return jumpFactor * func2->Eval(t, u)*n2[0];
		};
		// [[f2]]_y
		function<double(double, double)> jumpFunc2Y = [func2, n2, jumpFactor](double t, double u) {
			return jumpFactor * func2->Eval(t, u)*n2[1];
		};*/

		/*function<double(double, double)> functionToIntegrate = [meanFactor, meanGradFunc1X, meanGradFunc1Y, meanGradFunc2X, meanGradFunc2Y, jumpFunc1X, jumpFunc1Y, jumpFunc2X, jumpFunc2Y](double t, double u) {
			double meanGradFunc1_scal_jumpFunc2 = meanGradFunc1X(t, u)*jumpFunc2X(t, u) + meanGradFunc1Y(t, u)*jumpFunc2Y(t, u);
			double meanGradFunc2_scal_jumpFunc1 = meanGradFunc2X(t, u)*jumpFunc1X(t, u) + meanGradFunc2Y(t, u)*jumpFunc1Y(t, u);
			return meanGradFunc1_scal_jumpFunc2 + meanGradFunc2_scal_jumpFunc1;
		};*/

		function<double(double, double)> functionToIntegrate = [meanFactor, n1, n2, func1, func2, h](double t, double u) {
			double meanGradFunc1_scal_jumpFunc2 = meanFactor * (func1->EvalGradX(t, u) * n2[0] + func1->EvalGradY(t, u) * n2[1]) * func2->Eval(t, u);
			double meanGradFunc2_scal_jumpFunc1 = meanFactor * (func2->EvalGradX(t, u) * n1[0] + func2->EvalGradY(t, u) * n1[1]) * func1->Eval(t, u);
			return 2/h*(meanGradFunc1_scal_jumpFunc2 + meanGradFunc2_scal_jumpFunc1);
		};

		Element2DInterface* interf = (Element2DInterface*)interface;
		IPolynomialFunction* poly1 = dynamic_cast<IPolynomialFunction*>(func1);
		IPolynomialFunction* poly2 = dynamic_cast<IPolynomialFunction*>(func2);

		Square* square = (Square*)element1;

		GaussLegendre* gs = new GaussLegendre(poly1->GetDegree() + poly2->GetDegree() + 1);
		if (interf->X1 == interf->X2)
		{
			double t = interface == square->EastInterface ? 1 : -1;

			std::function<double(double)> func1D = [functionToIntegrate, t](double u) {
				return functionToIntegrate(t, u);
			};
			return -h / 2 * gs->Quadrature(func1D);
		}
		else if (interf->Y1 == interf->Y2)
		{
			double u = interface == square->NorthInterface ? 1 : -1;

			std::function<double(double)> func1D = [functionToIntegrate, u](double t) {
				return functionToIntegrate(t, u);
			};
			return -h / 2 * gs->Quadrature(func1D);
		}
		return 0;

		//double res = -2 / h * gs->Quadrature(functionToIntegrate, interf->X1, interf->X2, interf->Y1, interf->Y2);
		//return res;
		//return -interf->Integrate(functionToIntegrate);
	}

	double PenalizationTerm(ElementInterface* interface, Element* element1, IBasisFunction2D* func1, Element* element2, IBasisFunction2D* func2)
	{
		Square* square1 = static_cast<Square*>(element1);
		Square* square2 = static_cast<Square*>(element2);
		return this->PenalizationTerm(interface, square1, func1, square2, func2);
	}

	double PenalizationTerm(ElementInterface* interface, Square* element1, IBasisFunction2D* func1, Square* element2, IBasisFunction2D* func2)
	{
		//if (!interface->IsBetween(element1, element2))
		//	return 0;
		double h = element1->Width;
		auto n1 = element1->OuterNormalVector(interface);
		auto n2 = element2->OuterNormalVector(interface);

		//auto n2 = n1;
		/*double jumpFactor = interface->IsDomainBoundary ? 1 : -1;

		// [[f1]]_x
		function<double(double, double)> jumpFunc1X = [func1, n1, jumpFactor](double t, double u) {
			return jumpFactor * func1->Eval(t, u)*n1[0];
		};
		// [[f1]]_y
		function<double(double, double)> jumpFunc1Y = [func1, n1, jumpFactor](double t, double u) {
			return jumpFactor * func1->Eval(t, u)*n1[1];
		};
		// [[f2]]_x
		function<double(double, double)> jumpFunc2X = [func2, n2, jumpFactor](double t, double u) {
			return jumpFactor * func2->Eval(t, u)*n2[0];
		};
		// [[f2]]_y
		function<double(double, double)> jumpFunc2Y = [func2, n2, jumpFactor](double t, double u) {
			return jumpFactor * func2->Eval(t, u)*n2[1];
		};

		function<double(double, double)> functionToIntegrate = [jumpFunc1X, jumpFunc1Y, jumpFunc2X, jumpFunc2Y, h](double t, double u) {
			return jumpFunc1X(t, u)*jumpFunc2X(t, u) + jumpFunc1Y(t, u)*jumpFunc2Y(t, u);
		};*/

		function<double(double, double)> functionToIntegrate = [func1, func2, n1, n2, h](double t, double u) {
			return (n1[0] * n2[0] + n1[1] * n2[1]) * func1->Eval(t, u) * func2->Eval(t, u);
		};

		double integralJump1ScalarJump2 = 0;

		Element2DInterface* interf = (Element2DInterface*)interface;
		IPolynomialFunction* poly1 = dynamic_cast<IPolynomialFunction*>(func1);
		IPolynomialFunction* poly2 = dynamic_cast<IPolynomialFunction*>(func2);

		Square* square = (Square*)element1;

		GaussLegendre* gs = new GaussLegendre(poly1->GetDegree() + poly2->GetDegree() + 2);
		if (interf->X1 == interf->X2)
		{
			double t = interface == square->EastInterface ? 1 : -1;
			std::function<double(double)> func1D = [functionToIntegrate, t](double u) {
				return functionToIntegrate(t, u);
			};
			integralJump1ScalarJump2 = h / 2 * gs->Quadrature(func1D);
		}
		else if (interf->Y1 == interf->Y2)
		{
			double u = interface == square->NorthInterface ? 1 : -1;
			std::function<double(double)> func1D = [functionToIntegrate, u](double t) {
				return functionToIntegrate(t, u);
			};
			integralJump1ScalarJump2 = h / 2 * gs->Quadrature(func1D);
		}
		else
			return 0;
		//double integralJump1ScalarJump2 = /*h / 2 **/ gs->Quadrature(functionToIntegrate, interf->X1, interf->X2, interf->Y1, interf->Y2);
		//double integralJump1ScalarJump2 = interf->Integrate(functionToIntegrate);

		return this->_penalizationCoefficient * integralJump1ScalarJump2;
	}

	double RightHandSide(Element* element, IBasisFunction2D* func)
	{
		Square* square = static_cast<Square*>(element);
		return this->RightHandSide(square, func);
	}

	double RightHandSide(Square* element, IBasisFunction2D* func)
	{
		double x1 = element->X;
		double x2 = element->X + element->Width;
		double y1 = element->Y;
		double y2 = element->Y + element->Width;

		function<double(double, double)> sourceTimesBasisFunction = [this, func, x1, x2, y1, y2](double t, double u) {
			return this->_sourceFunction((x2 - x1) / 2 * t + (x2 + x1) / 2, (y2 - y1) / 2 * u + (y2 + y1) / 2) * func->Eval(t, u);
		};

		GaussLegendre* gs = new GaussLegendre();
		return (x2 - x1) * (y2 - y1) / 4 * gs->Quadrature(sourceTimesBasisFunction);
	}
};