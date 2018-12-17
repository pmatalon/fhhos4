#pragma once
#include <map>
#include <functional>
#include <math.h>
#include "FunctionalBasisWithNumbers.h"
#include "CartesianGrid1D.h"
#include "IBasisFunction1D.h"
//#include "Utils.h"
using namespace std;

class FunctionalGlobalBasis1D : public FunctionalBasisWithNumbers
{
protected:
	CartesianGrid1D* _grid;
	int _penalizationCoefficient;
	function<double(double)> _sourceFunction;

	map<int, IBasisFunction1D*> _localFunctions;

public:
	FunctionalGlobalBasis1D(CartesianGrid1D* grid, int penalizationCoefficient, function<double(double)> sourceFunction)
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
		function<double(double)> functionToIntegrate = [func1, func2](double x) {
			return func1->EvalDerivative(x)*func2->EvalDerivative(x);
		};

		GaussLegendre* gs = new GaussLegendre(func1->GetDegree() + func2->GetDegree());
		return gs->Quadrature(functionToIntegrate, this->_grid->XLeft(element), this->_grid->XRight(element));
		//return Utils::Integral(functionToIntegrate, this->_grid->XLeft(element), this->_grid->XRight(element));
	}

	double CouplingTerm(BigNumber interface, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2)
	{
		if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;

		return MeanGrad(element1, func1, interface) * Jump(element2, func2, interface) + MeanGrad(element2, func2, interface) * Jump(element1, func1, interface);
	}

	double PenalizationTerm(BigNumber point, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2)
	{
		if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;

		return this->_penalizationCoefficient * Jump(element1, func1, point) * Jump(element2, func2, point);
	}

	double RightHandSide(BigNumber element, IBasisFunction1D* func)
	{
		function<double(double)> sourceTimesBasisFunction = [this, func](double x) {
			return this->_sourceFunction(x) * func->Eval(x);
		};

		GaussLegendre* gs = new GaussLegendre();
		return gs->Quadrature(sourceTimesBasisFunction, this->_grid->XLeft(element), this->_grid->XRight(element));
		//return Utils::Integral(sourceTimesBasisFunction, this->_grid->XLeft(element), this->_grid->XRight(element));
	}

	double MeanGrad(BigNumber element, IBasisFunction1D* func, BigNumber interface)
	{
		if (this->_grid->IsBoundaryLeft(interface) || this->_grid->IsBoundaryRight(interface))
			return func->EvalDerivative(this->_grid->X(interface));
		return 0.5 * (func->EvalDerivative(this->_grid->X(interface)));
	}

	double Jump(BigNumber element, IBasisFunction1D* func, BigNumber interface)
	{
		int factor = this->_grid->IsLeftInterface(element, interface) ? 1 : -1;
		return factor * (func->Eval(this->_grid->X(interface)));
	}
};