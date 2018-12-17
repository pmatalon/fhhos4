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
#include "IPolynomialFunction.h"
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
			return func1->EvalGradX(t, u)*func2->EvalGradX(t, u) + func1->EvalGradY(t, u)*func2->EvalGradY(t, u);
		};

		IPolynomialFunction* poly1 = dynamic_cast<IPolynomialFunction*>(func1);
		IPolynomialFunction* poly2 = dynamic_cast<IPolynomialFunction*>(func2);

		GaussLegendre gs(poly1->GetDegree() + poly2->GetDegree());
		//return 4 / pow(h, 2) * gs.Quadrature(functionToIntegrate);
		//return pow(h, 2) / 4 * gs.Quadrature(functionToIntegrate);
		return gs.Quadrature(functionToIntegrate);
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

		Element2DInterface* interf = (Element2DInterface*)interface;
		IPolynomialFunction* poly1 = dynamic_cast<IPolynomialFunction*>(func1);
		IPolynomialFunction* poly2 = dynamic_cast<IPolynomialFunction*>(func2);

		function<double(double, double)> functionToIntegrate = [meanFactor, n1, n2, func1, func2, h, element1, element2, poly1, poly2](double t, double u) {
			double meanGradFunc1_scal_jumpFunc2 = meanFactor * (func1->EvalGradX(t, u) * n2[0] + func1->EvalGradY(t, u) * n2[1]) * func2->Eval(t, u);
			double meanGradFunc2_scal_jumpFunc1 = meanFactor * (func2->EvalGradX(t, u) * n1[0] + func2->EvalGradY(t, u) * n1[1]) * func1->Eval(t, u);
			/*if (element1->Number == 1 && element2->Number == 1 && poly1->GetDegree() == 1 && poly2->GetDegree() == 1)
			{
				cout << func2->Eval(t, u) << endl;
			}*/
			return 2/h*(meanGradFunc1_scal_jumpFunc2 + meanGradFunc2_scal_jumpFunc1);
		};


		Square* square = (Square*)element1;

		GaussLegendre gs(poly1->GetDegree() + poly2->GetDegree() + 1);
		std::function<double(double)> func1D;
		if (interf->IsVertical())
		{
			double t = interface == square->EastInterface ? 1 : -1;
			func1D = [functionToIntegrate, t](double u) {
				return functionToIntegrate(t, u);
			};
		}
		else if (interf->IsHorizontal())
		{
			double u = interface == square->NorthInterface ? 1 : -1;
			func1D = [functionToIntegrate, u](double t) {
				return functionToIntegrate(t, u);
			};
		}
		else
			return 0;

		return -h / 2 * gs.Quadrature(func1D);
		//return -pow(h, 2) / 4 * gs.Quadrature(func1D);
		//return -2/h* gs.Quadrature(func1D);
		//return -gs.Quadrature(func1D);
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
		
		function<double(double, double)> functionToIntegrate = [func1, func2, n1, n2, h](double t, double u) {
			return (n1[0] * n2[0] + n1[1] * n2[1]) * func1->Eval(t, u) * func2->Eval(t, u);
		};

		double integralJump1ScalarJump2 = 0;

		Element2DInterface* interf = (Element2DInterface*)interface;
		IPolynomialFunction* poly1 = dynamic_cast<IPolynomialFunction*>(func1);
		IPolynomialFunction* poly2 = dynamic_cast<IPolynomialFunction*>(func2);

		Square* square = (Square*)element1;

		GaussLegendre gs(poly1->GetDegree() + poly2->GetDegree() + 2);
		std::function<double(double)> func1D;
		if (interf->IsVertical())
		{
			double t = interface == square->EastInterface ? 1 : -1;
			func1D = [functionToIntegrate, t](double u) {
				return functionToIntegrate(t, u);
			};
		}
		else if (interf->IsHorizontal())
		{
			double u = interface == square->NorthInterface ? 1 : -1;
			func1D = [functionToIntegrate, u](double t) {
				return functionToIntegrate(t, u);
			};
		}
		else
			return 0;

		integralJump1ScalarJump2 = h / 2 * gs.Quadrature(func1D);
		//integralJump1ScalarJump2 = pow(h, 2) / 4 * gs.Quadrature(func1D);
		//integralJump1ScalarJump2 = 2 / h * gs.Quadrature(func1D);
		//integralJump1ScalarJump2 = gs.Quadrature(func1D);
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

		GaussLegendre gs;
		return (x2 - x1) * (y2 - y1) / 4 * gs.Quadrature(sourceTimesBasisFunction);
	}
};