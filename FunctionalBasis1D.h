#pragma once
#include <map>
#include <functional>
#include <math.h>
#include "FunctionalBasisWithNumbers.h"
#include "CartesianGrid1D.h"
#include "BasisFunction1D.h"
#include "Utils.h"
using namespace std;

class FunctionalBasis1D : public FunctionalBasisWithNumbers
{
protected:
	CartesianGrid1D* _grid;
	int _penalizationCoefficient;
	function<double(double)> _sourceFunction;

	map<int, BasisFunction1D*> _localFunctions;

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

	BasisFunction1D* GetLocalBasisFunction(BigNumber element, int localFunctionNumber)
	{
		return this->_localFunctions[localFunctionNumber];
	}

	BigNumber GlobalFunctionNumber(BigNumber element, int localFunctionNumber)
	{
		return element * NumberOfLocalFunctionsInElement(0) + localFunctionNumber + 1; // +1 so that the numbers start at 1
	}

	double VolumicTerm(BigNumber element, BasisFunction1D* func1, BasisFunction1D* func2)
	{
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		function<double(double)> functionToIntegrate = [func1, func2](double x) {
			return func1->EvalDerivative(x)*func2->EvalDerivative(x);
		};

		return 2/(b - a) * Utils::IntegralOnReferenceInterval(functionToIntegrate);
	}

	double CouplingTerm(BigNumber interface, BigNumber element1, BasisFunction1D* func1, BigNumber element2, BasisFunction1D* func2)
	{
		if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;

		return MeanGrad(element1, func1, interface) * Jump(element2, func2, interface) + MeanGrad(element2, func2, interface) * Jump(element1, func1, interface);
	}

	double PenalizationTerm(BigNumber point, BigNumber element1, BasisFunction1D* func1, BigNumber element2, BasisFunction1D* func2)
	{
		if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;

		return this->_penalizationCoefficient * Jump(element1, func1, point) * Jump(element2, func2, point);
	}

	double RightHandSide(BigNumber element, BasisFunction1D* func)
	{
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		function<double(double)> sourceTimesBasisFunction = [this, func, a, b](double x) {
			return this->_sourceFunction((b-a)/2*x + (a+b)/2) * func->Eval(x);
		};
		return (b - a) / 2 * Utils::IntegralOnReferenceInterval(sourceTimesBasisFunction);
	}

	double MeanGrad(BigNumber element, BasisFunction1D* func, BigNumber interface)
	{
		int xRef = this->_grid->IsLeftInterface(element, interface) ? -1 : 1; // xRef in [-1, 1]
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		if (this->_grid->IsBoundaryLeft(interface) || this->_grid->IsBoundaryRight(interface))
			return 2 / (b - a) * func->EvalDerivative(xRef);
		return 1 / (b - a) * (func->EvalDerivative(xRef));
	}

	double Jump(BigNumber element, BasisFunction1D* func, BigNumber interface)
	{
		int xRef = this->_grid->IsLeftInterface(element, interface) ? -1 : 1; // xRef in [-1, 1]
		int factor = this->_grid->IsLeftInterface(element, interface) ? 1 : -1;
		return factor * (func->Eval(xRef));
	}
};