#pragma once
#include <functional>
#include <math.h>
#include "IPoisson1D_DGTerms_OLD.h"
#include "Utils.h"
#include "Element.h"
#include "CartesianGrid1DOLD.h"
using namespace std;

class Poisson1D_DGTerms_LocalBasisOLD : public IPoisson1D_DGTerms
{
protected:
	CartesianGrid1DOLD* _grid;
	function<double(double)> _sourceFunction;

public:
	Poisson1D_DGTerms_LocalBasisOLD(CartesianGrid1DOLD* grid, function<double(double)> sourceFunction)
	{
		this->_grid = grid;
		this->_sourceFunction = sourceFunction;
	}

	bool IsGlobalBasis() { return false; }

	double VolumicTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2)
	{
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		int nQuadPoints = func1->GetDegree() + func2->GetDegree();
		// defined on [-1, 1]
		function<double(double)> functionToIntegrate = [func1, func2](double t) {
			return func1->EvalDerivative(t)*func2->EvalDerivative(t);
		};

		GaussLegendre gs(nQuadPoints);
		return 2 / (b - a) * gs.Quadrature(functionToIntegrate);
	}

	double MassTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2)
	{
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		GaussLegendre gs(func1->GetDegree() + func2->GetDegree());

		// defined on [-1, 1]
		function<double(double)> functionToIntegrate = [func1, func2](double t) {
			return func1->Eval(t)*func2->Eval(t);
		};

		return (b - a) / 2 * gs.Quadrature(functionToIntegrate);
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

		GaussLegendre gs;

		function<double(double)> sourceTimesBasisFunction = [this, func, a, b](double t) {
			return this->_sourceFunction((b - a) / 2 * t + (a + b) / 2) * func->Eval(t);
		};

		return (b - a) / 2 * gs.Quadrature(sourceTimesBasisFunction);
	}

	double MeanDerivative(BigNumber element, IBasisFunction1D* func, BigNumber interface)
	{
		double t = this->_grid->IsLeftInterface(element, interface) ? -1 : 1; // t in [-1, 1]
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		double jacobian = 2 / (b - a);

		if (this->_grid->IsBoundaryLeft(interface) || this->_grid->IsBoundaryRight(interface))
			return jacobian * func->EvalDerivative(t);
		return 0.5 * jacobian * func->EvalDerivative(t);
	}

	double Jump(BigNumber element, IBasisFunction1D* func, BigNumber interface)
	{
		double t = this->_grid->IsLeftInterface(element, interface) ? -1 : 1; // t in [-1, 1]
		int factor = this->_grid->IsLeftInterface(element, interface) ? 1 : -1;
		return factor * (func->Eval(t));
	}
};