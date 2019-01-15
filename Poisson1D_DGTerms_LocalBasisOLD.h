#pragma once
#include <functional>
#include <math.h>
#include "IPoisson1D_DGTerms.h"
#include "Utils.h"
#include "Element.h"
#include "CartesianGrid1D.h"
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
		if (func1->DefinitionInterval().Left == -1 && func1->DefinitionInterval().Right == 1)
		{
			// defined on [-1, 1]
			function<double(double)> functionToIntegrate = [func1, func2](double t) {
				return func1->EvalDerivative(t)*func2->EvalDerivative(t);
			};

			GaussLegendre gs(nQuadPoints);
			return 2 / (b - a) * gs.Quadrature(functionToIntegrate);
		}
		else
		{
			function<double(double)> functionToIntegrate = [func1, func2](double u) {
				//double u = 0.5 * t + 0.5;
				return func1->EvalDerivative(u)*func2->EvalDerivative(u);
			};

			return 1 / (b - a) * Utils::Integral(nQuadPoints, functionToIntegrate, 0, 1);
		}
	}

	double MassTerm(BigNumber element, IBasisFunction1D* func1, IBasisFunction1D* func2)
	{
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		GaussLegendre gs(func1->GetDegree() + func2->GetDegree());

		if (func1->DefinitionInterval().Left == -1 && func1->DefinitionInterval().Right == 1)
		{
			// defined on [-1, 1]
			function<double(double)> functionToIntegrate = [func1, func2](double t) {
				return func1->Eval(t)*func2->Eval(t);
			};

			return (b - a) / 2 * gs.Quadrature(functionToIntegrate);
		}
		else
		{
			function<double(double)> functionToIntegrate = [func1, func2](double u) {
				//double u = 0.5 * t + 0.5;
				return func1->Eval(u)*func2->Eval(u);
			};

			return (b - a) * Utils::Integral(func1->GetDegree() + func2->GetDegree(), functionToIntegrate, 0, 1);
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

		if (func->DefinitionInterval().Left == -1 && func->DefinitionInterval().Right == 1)
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
		double t = this->_grid->IsLeftInterface(element, interface) ? func->DefinitionInterval().Left : func->DefinitionInterval().Right; // t in [-1, 1]
		double a = this->_grid->XLeft(element);
		double b = this->_grid->XRight(element);

		double jacobian = 0;
		if (func->DefinitionInterval().Left == -1 && func->DefinitionInterval().Right == 1)
			jacobian = 2 / (b - a);
		else
			jacobian = 1 / (b - a);

		if (this->_grid->IsBoundaryLeft(interface) || this->_grid->IsBoundaryRight(interface))
			return jacobian * func->EvalDerivative(t);
		return 0.5 * jacobian * func->EvalDerivative(t);
	}

	double Jump(BigNumber element, IBasisFunction1D* func, BigNumber interface)
	{
		double t = this->_grid->IsLeftInterface(element, interface) ? func->DefinitionInterval().Left : func->DefinitionInterval().Right; // t in [-1, 1]
		int factor = this->_grid->IsLeftInterface(element, interface) ? 1 : -1;
		return factor * (func->Eval(t));
	}
};