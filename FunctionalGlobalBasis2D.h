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
using namespace std;

class FunctionalGlobalBasis2D : public FunctionalBasisWithObjects
{
protected:
	CartesianGrid2D* _grid;
	int _penalizationCoefficient;
	function<double(double, double)> _sourceFunction;

	map<int, IBasisFunction2D*> _localFunctions;

public:
	FunctionalGlobalBasis2D(CartesianGrid2D* grid, int penalizationCoefficient, function<double(double, double)> sourceFunction)
	{
		this->_grid = grid;
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
		function<double(double, double)> functionToIntegrate = [func1, func2](double x, double y) {
			return func1->EvalGradX(x, y)*func2->EvalGradX(x, y) + func1->EvalGradY(x, y)*func2->EvalGradY(x, y);
		};

		return Utils::Integral(functionToIntegrate, element->X, element->X + element->Width, element->Y, element->Y + element->Width);
	}

	double CouplingTerm(ElementInterface* interface, Element* element1, IBasisFunction2D* func1, Element* element2, IBasisFunction2D* func2)
	{
		if (!interface->IsBetween(element1, element2))
			return 0;

		auto n1 = element1->OuterNormalVector(interface);
		auto n2 = element2->OuterNormalVector(interface);

		double meanFactor = interface->IsDomainBoundary ? 1 : 0.5;
		double jumpFactor = interface->IsDomainBoundary ? 1 : -1;

		// {{grad f1}}_x
		function<double(double, double)> meanGradFunc1X = [func1, meanFactor](double x, double y) {
			return meanFactor * func1->EvalGradX(x, y);
		};
		// {{grad f1}}_y
		function<double(double, double)> meanGradFunc1Y = [func1, meanFactor](double x, double y) {
			return meanFactor * func1->EvalGradY(x, y);
		};
		// {{grad f2}}_x
		function<double(double, double)> meanGradFunc2X = [func2, meanFactor](double x, double y) {
			return meanFactor * func2->EvalGradX(x, y);
		};
		// {{grad f2}}_y
		function<double(double, double)> meanGradFunc2Y = [func2, meanFactor](double x, double y) {
			return meanFactor * func2->EvalGradY(x, y);
		};
		// [[f1]]_x
		function<double(double, double)> jumpFunc1X = [func1, n1, jumpFactor](double x, double y) {
			return jumpFactor * func1->Eval(x, y)*n1[0];
		};
		// [[f1]]_y
		function<double(double, double)> jumpFunc1Y = [func1, n1, jumpFactor](double x, double y) {
			return jumpFactor * func1->Eval(x, y)*n1[1];
		};
		// [[f2]]_x
		function<double(double, double)> jumpFunc2X = [func2, n2, jumpFactor](double x, double y) {
			return jumpFactor * func2->Eval(x, y)*n2[0];
		};
		// [[f2]]_y
		function<double(double, double)> jumpFunc2Y = [func2, n2, jumpFactor](double x, double y) {
			return jumpFactor * func2->Eval(x, y)*n2[1];
		};

		function<double(double, double)> functionToIntegrate = [meanGradFunc1X, meanGradFunc1Y, meanGradFunc2X, meanGradFunc2Y, jumpFunc1X, jumpFunc1Y, jumpFunc2X, jumpFunc2Y](double x, double y) {
			double meanGradF1X = meanGradFunc1X(x, y);
			double meanGradF1Y = meanGradFunc1Y(x, y);
			double meanGradF2X = meanGradFunc2X(x, y);
			double meanGradF2Y = meanGradFunc2Y(x, y);
			double jumpF1X = jumpFunc1X(x, y);
			double jumpF1Y = jumpFunc1Y(x, y);
			double jumpF2X = jumpFunc2X(x, y);
			double jumpF2Y = jumpFunc2Y(x, y);

			double meanGradFunc1_scal_jumpFunc2 = meanGradF1X*jumpF2X + meanGradF1Y*jumpF2Y;
			double meanGradFunc2_scal_jumpFunc1 = meanGradF2X*jumpF1X + meanGradF2Y*jumpF1Y;
			/*double meanGradFunc1_scal_jumpFunc2 = meanGradFunc1X(x, y)*jumpFunc2X(x, y) + meanGradFunc1Y(x, y)*jumpFunc2Y(x, y);
			double meanGradFunc2_scal_jumpFunc1 = meanGradFunc2X(x, y)*jumpFunc1X(x, y) + meanGradFunc2Y(x, y)*jumpFunc1Y(x, y);*/
			return meanGradFunc1_scal_jumpFunc2 + meanGradFunc2_scal_jumpFunc1;
		};

		Element2DInterface* interf = (Element2DInterface*)interface;
		return -interf->Integrate(functionToIntegrate);
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

		auto n1 = element1->OuterNormalVector(interface);
		auto n2 = element2->OuterNormalVector(interface);
		double jumpFactor = interface->IsDomainBoundary ? 1 : -1;

		// [[f1]]_x
		function<double(double, double)> jumpFunc1X = [func1, n1, jumpFactor](double x, double y) {
			return jumpFactor * func1->Eval(x, y)*n1[0];
		};
		// [[f1]]_y
		function<double(double, double)> jumpFunc1Y = [func1, n1, jumpFactor](double x, double y) {
			return jumpFactor * func1->Eval(x, y)*n1[1];
		};
		// [[f2]]_x
		function<double(double, double)> jumpFunc2X = [func2, n2, jumpFactor](double x, double y) {
			return jumpFactor * func2->Eval(x, y)*n2[0];
		};
		// [[f2]]_y
		function<double(double, double)> jumpFunc2Y = [func2, n2, jumpFactor](double x, double y) {
			return jumpFactor * func2->Eval(x, y)*n2[1];
		};

		function<double(double, double)> functionToIntegrate = [jumpFunc1X, jumpFunc1Y, jumpFunc2X, jumpFunc2Y](double x, double y) {
			return jumpFunc1X(x, y)*jumpFunc2X(x, y) + jumpFunc1Y(x, y)*jumpFunc2Y(x, y);
		};
		
		Element2DInterface* interf = (Element2DInterface*)interface;
		double integralJump1ScalarJump2 = interf->Integrate(functionToIntegrate);

		return this->_penalizationCoefficient * integralJump1ScalarJump2;
	}

	double RightHandSide(Element* element, IBasisFunction2D* func)
	{
		Square* square = static_cast<Square*>(element);
		return this->RightHandSide(square, func);
	}

	double RightHandSide(Square* element, IBasisFunction2D* func)
	{
		function<double(double, double)> sourceTimesBasisFunction = [this, func](double x, double y) {
			return this->_sourceFunction(x, y) * func->Eval(x, y);
		};
		return Utils::Integral(sourceTimesBasisFunction, element->X, element->X + element->Width, element->Y, element->Y + element->Width);
	}
};