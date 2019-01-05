#pragma once
#include <functional>
#include <math.h>
#include "IPoisson1D_DGTerms.h"
#include "Utils.h"
#include "Element.h"
#include "CartesianGrid1D.h"
using namespace std;

class Poisson1D_DGTerms_GlobalBasis : public IPoisson1D_DGTerms
{
protected:
	CartesianGrid1D* _grid;
	function<double(double)> _sourceFunction;

public:
	Poisson1D_DGTerms_GlobalBasis(CartesianGrid1D* grid, function<double(double)> sourceFunction)
	{
		this->_grid = grid;
		this->_sourceFunction = sourceFunction;
	}

	bool IsGlobalBasis() { return true; }

	double VolumicTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2)
	{
		function<double(double)> functionToIntegrate = [func1, func2](double x) {
			return func1->EvalDerivative(x)*func2->EvalDerivative(x);
		};

		GaussLegendre* gs = new GaussLegendre(func1->GetDegree() + func2->GetDegree());
		return gs->Quadrature(functionToIntegrate, this->_grid->XLeft(element), this->_grid->XRight(element));
		//return Utils::Integral(functionToIntegrate, this->_grid->XLeft(element), this->_grid->XRight(element));
	}

	double MassTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2)
	{
		function<double(double)> functionToIntegrate = [func1, func2](double x) {
			return func1->Eval(x)*func2->Eval(x);
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

	double PenalizationTerm(BigNumber point, BigNumber element1, IBasisFunction1D* func1, BigNumber element2, IBasisFunction1D* func2, double penalizationCoefficient)
	{
		if (element2 > element1 + 1 || element1 > element2 + 1)
			return 0;

		return penalizationCoefficient * Jump(element1, func1, point) * Jump(element2, func2, point);
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