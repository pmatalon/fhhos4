#pragma once
#include <map>
#include <functional>
#include <math.h>
#include "FunctionalBasisWithNumbers.h"
#include "CartesianGrid1D.h"
#include "IBasisFunction1D.h"
//#include "Utils.h"
#include "GaussLegendre.h"
#include "IPolynomialFunction.h"
using namespace std;

class FunctionalBasis1D : public FunctionalBasisWithNumbers
{
protected:
	CartesianGrid1D* _grid;
	int _penalizationCoefficient;
	function<double(double)> _sourceFunction;

	map<int, IBasisFunction1D*> _localFunctions;

public:
	FunctionalBasis1D(CartesianGrid1D* grid, int penalizationCoefficient, function<double(double)> sourceFunction)
	{
		this->_grid = grid;
		this->_penalizationCoefficient = penalizationCoefficient;
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

		function<double(double)> functionToIntegrate = [func1, func2](double t) {
			return func1->EvalDerivative(t)*func2->EvalDerivative(t);
		};

		IPolynomialFunction* poly1 = dynamic_cast<IPolynomialFunction*>(func1);
		IPolynomialFunction* poly2 = dynamic_cast<IPolynomialFunction*>(func2);

		GaussLegendre* gs = new GaussLegendre(poly1->GetDegree() + poly2->GetDegree());
		return 2 / (b - a) * gs->Quadrature(functionToIntegrate);
		//return 2/(b - a) * Utils::IntegralOnReferenceInterval(functionToIntegrate);
	}

	double CouplingTerm(BigNumber interface, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2)
	{
		if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;

		return MeanDerivative(element1, func1, interface) * Jump(element2, func2, interface) + MeanDerivative(element2, func2, interface) * Jump(element1, func1, interface);
	}

	double PenalizationTerm(BigNumber point, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2)
	{
		if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;

		return this->_penalizationCoefficient * Jump(element1, func1, point) * Jump(element2, func2, point);
	}

	double RightHandSide(BigNumber element, IBasisFunction1D* func)
	{
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		function<double(double)> sourceTimesBasisFunction = [this, func, a, b](double t) {
			return this->_sourceFunction((b-a)/2*t + (a+b)/2) * func->Eval(t);
		};

		GaussLegendre* gs = new GaussLegendre();
		return (b - a) / 2 * gs->Quadrature(sourceTimesBasisFunction);
		//return (b - a) / 2 * Utils::IntegralOnReferenceInterval(sourceTimesBasisFunction);
	}

	double MeanDerivative(BigNumber element, IBasisFunction1D* func, BigNumber interface)
	{
		int t = this->_grid->IsLeftInterface(element, interface) ? -1 : 1; // t in [-1, 1]
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		if (this->_grid->IsBoundaryLeft(interface) || this->_grid->IsBoundaryRight(interface))
			return 2 / (b - a) * func->EvalDerivative(t);
		return 1 / (b - a) * (func->EvalDerivative(t));
	}

	double Jump(BigNumber element, IBasisFunction1D* func, BigNumber interface)
	{
		int t = this->_grid->IsLeftInterface(element, interface) ? -1 : 1; // t in [-1, 1]
		int factor = this->_grid->IsLeftInterface(element, interface) ? 1 : -1;
		return factor * (func->Eval(t));
	}
};