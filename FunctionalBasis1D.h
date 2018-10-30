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

	BigNumber GlobalFunctionNumber(BigNumber element, int localFunctionNumber)
	{
		return element * NumberOfLocalFunctionsInElement(0) + localFunctionNumber + 1; // +1 so that the numbers start at 1
	}

	//virtual double VolumicTerm(BigNumber element, int localFunctionNumber1, int localFunctionNumber2) = 0;
	double VolumicTerm(BigNumber element, int localFunctionNumber1, int localFunctionNumber2)
	{
		BasisFunction1D* func1 = this->_localFunctions[localFunctionNumber1];
		BasisFunction1D* func2 = this->_localFunctions[localFunctionNumber2];

		function<double(double)> functionToIntegrate = [func1, func2](double x) {
			return func1->EvalGrad(x)*func2->EvalGrad(x);
		};

		return Utils::Integral(functionToIntegrate, this->_grid->XLeft(element), this->_grid->XRight(element));
	}

	double CouplingTerm(BigNumber interface, BigNumber element1, int localFunctionNumber1, BigNumber element2, int localFunctionNumber2)
	{
		if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;

		BasisFunction1D* func1 = this->_localFunctions[localFunctionNumber1];
		BasisFunction1D* func2 = this->_localFunctions[localFunctionNumber2];

		return MeanGrad(element1, func1, interface) * Jump(element2, func2, interface) + MeanGrad(element2, func2, interface) * Jump(element1, func1, interface);
	}

	double PenalizationTerm(BigNumber point, BigNumber element1, int localFunctionNumber1, BigNumber element2, int localFunctionNumber2)
	{
		if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;

		BasisFunction1D* func1 = this->_localFunctions[localFunctionNumber1];
		BasisFunction1D* func2 = this->_localFunctions[localFunctionNumber2];

		return this->_penalizationCoefficient * Jump(element1, func1, point) * Jump(element2, func2, point);
	}

	double RightHandSide(BigNumber element, int localFunctionNumber)
	{
		BasisFunction1D* func1 = this->_localFunctions[localFunctionNumber];
		function<double(double)> sourceTimesBasisFunction = [this, func1](double x) {
			return this->_sourceFunction(x) * func1->Eval(x);
		};
		return Utils::Integral(sourceTimesBasisFunction, this->_grid->XLeft(element), this->_grid->XRight(element));
	}

	double MeanGrad(BigNumber element, BasisFunction1D* func, BigNumber interface)
	{
		if (this->_grid->IsBoundaryLeft(interface) || this->_grid->IsBoundaryRight(interface))
			return func->EvalGrad(this->_grid->X(interface));
		return 0.5 * (func->EvalGrad(this->_grid->X(interface)));
	}

	double Jump(BigNumber element, BasisFunction1D* func, BigNumber interface)
	{
		int factor = this->_grid->IsLeftInterface(element, interface) ? 1 : -1;
		return factor * (func->Eval(this->_grid->X(interface)));
	}
};