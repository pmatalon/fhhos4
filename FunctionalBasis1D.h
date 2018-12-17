#pragma once
#include <map>
#include <functional>
#include <math.h>
#include "FunctionalBasisWithNumbers.h"
#include "CartesianGrid1D.h"
#include "IBasisFunction1D.h"
#include "GaussLegendre.h"
#include "IPolynomialFunction.h"
#include "Utils.h"
using namespace std;

class FunctionalBasis1D : public FunctionalBasisWithNumbers
{
protected:
	CartesianGrid1D* _grid;
	function<double(double)> _sourceFunction;

	map<int, IBasisFunction1D*> _localFunctions;

public:
	FunctionalBasis1D(CartesianGrid1D* grid, function<double(double)> sourceFunction)
	{
		this->_grid = grid;
		this->_sourceFunction = sourceFunction;
	}

	int NumberOfLocalFunctionsInElement(BigNumber element)
	{
		return static_cast<int>(this->_localFunctions.size());
	}

	IBasisFunction1D* GetLocalBasisFunction(BigNumber element, int localFunctionNumber)
	{
		return this->_localFunctions[localFunctionNumber];
	}

	BigNumber GlobalFunctionNumber(BigNumber element, int localFunctionNumber)
	{
		return element * NumberOfLocalFunctionsInElement(0) + localFunctionNumber + 1; // +1 so that the numbers start at 1
	}

	double VolumicTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2)
	{
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		IPolynomialFunction* poly1 = dynamic_cast<IPolynomialFunction*>(func1);
		IPolynomialFunction* poly2 = dynamic_cast<IPolynomialFunction*>(func2);
		GaussLegendre gs(poly1->GetDegree() + poly2->GetDegree());

		if (func1->ReferenceInterval().Left == -1 && func1->ReferenceInterval().Right == 1)
		{
			// defined on [-1, 1]
			function<double(double)> functionToIntegrate = [func1, func2](double t) {
				return func1->EvalDerivative(t)*func2->EvalDerivative(t);
			};

			return 2 / (b - a) * gs.Quadrature(functionToIntegrate);
		}
		else
		{
			function<double(double)> functionToIntegrate = [func1, func2](double u) {
				//double u = 0.5 * t + 0.5;
				return func1->EvalDerivative(u)*func2->EvalDerivative(u);
			};

			return 1 / (b - a) * Utils::Integral(poly1->GetDegree() + poly2->GetDegree(), functionToIntegrate, 0, 1);
		}
	}

	double MassTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2)
	{
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		IPolynomialFunction* poly1 = dynamic_cast<IPolynomialFunction*>(func1);
		IPolynomialFunction* poly2 = dynamic_cast<IPolynomialFunction*>(func2);
		GaussLegendre gs(poly1->GetDegree() + poly2->GetDegree());

		if (func1->ReferenceInterval().Left == -1 && func1->ReferenceInterval().Right == 1)
		{
			// defined on [-1, 1]
			function<double(double)> functionToIntegrate = [func1, func2](double t) {
				return func1->Eval(t)*func2->Eval(t);
			};

			return 2 / (b - a) * gs.Quadrature(functionToIntegrate);
		}
		else
		{
			function<double(double)> functionToIntegrate = [func1, func2](double u) {
				//double u = 0.5 * t + 0.5;
				return func1->Eval(u)*func2->Eval(u);
			};

			return 1 / (b - a) * Utils::Integral(poly1->GetDegree() + poly2->GetDegree(), functionToIntegrate, 0, 1);
		}
	}

	double CouplingTerm(BigNumber interface, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2)
	{
		if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;

		return MeanDerivative(element1, func1, interface) * Jump(element2, func2, interface) + MeanDerivative(element2, func2, interface) * Jump(element1, func1, interface);
	}

	double PenalizationTerm(BigNumber point, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2, double penalizationCoefficient)
	{
		if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;

		return penalizationCoefficient * Jump(element1, func1, point) * Jump(element2, func2, point);
	}

	double RightHandSide(BigNumber element, IBasisFunction1D* func)
	{
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

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

	double MeanDerivative(BigNumber element, IBasisFunction1D* func, BigNumber interface)
	{
		double t = this->_grid->IsLeftInterface(element, interface) ? func->ReferenceInterval().Left : func->ReferenceInterval().Right; // t in [-1, 1]
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		double jacobian = 0;
		if (func->ReferenceInterval().Left == -1 && func->ReferenceInterval().Right == 1)
			jacobian = 2 / (b - a);
		else
			jacobian = 1 / (b - a);

		if (this->_grid->IsBoundaryLeft(interface) || this->_grid->IsBoundaryRight(interface))
			return jacobian * func->EvalDerivative(t);
		return 0.5 * jacobian * func->EvalDerivative(t);
	}

	double Jump(BigNumber element, IBasisFunction1D* func, BigNumber interface)
	{
		double t = this->_grid->IsLeftInterface(element, interface) ? func->ReferenceInterval().Left : func->ReferenceInterval().Right; // t in [-1, 1]
		int factor = this->_grid->IsLeftInterface(element, interface) ? 1 : -1;
		return factor * (func->Eval(t));
	}
};