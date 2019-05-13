#pragma once
#include "Element.h"
#include "CartesianShape.h"

template <int Dim>
class CartesianElement : virtual public Element<Dim>, public CartesianShape<Dim>
{
public:

	CartesianElement(BigNumber number, DomPoint origin, double width) : Element<Dim>(number), CartesianShape<Dim>(origin, width)
	{
	}

	double GetDiameter() override
	{
		return CartesianShape<Dim>::Width;
	}

	double Measure()
	{
		return CartesianShape<Dim>::Measure();
	}

	double DiffusionCoefficient(DiffusionPartition diffusionPartition) override
	{
		DomPoint origin = CartesianShape<Dim>::Origin;
		return diffusionPartition.Coefficient(origin);
	}

	double Integral(BasisFunction<Dim>* phi)
	{
		return CartesianShape<Dim>::Integral(phi);
	}

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return CartesianShape<Dim>::MassTerm(phi1, phi2);
	}

	double IntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return CartesianShape<Dim>::IntegralGradGrad(phi1, phi2);
	}

	double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return CartesianShape<Dim>::ComputeIntegralGradGrad(phi1, phi2);
	}

	DomPoint ConvertToDomain(RefPoint referenceElementPoint)
	{
		return CartesianShape<Dim>::ConvertToDomain(referenceElementPoint);
	}

	RefPoint ConvertToReference(DomPoint domainPoint)
	{
		return CartesianShape<Dim>::ConvertToReference(domainPoint);
	}

	vector<RefPoint> GetNodalPoints(FunctionalBasis<Dim>* basis)
	{
		return CartesianShape<Dim>::GetNodalPoints(basis);
	}

	double SourceTerm(BasisFunction<Dim>* phi, SourceFunction* f)
	{
		function<double(RefPoint)> sourceTimesBasisFunction = [this, f, phi](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return f->Eval(domainPoint) * phi->Eval(refElementPoint);
		};

		CartesianShape<Dim>::ComputeIntegral(sourceTimesBasisFunction, 0);
	}

	double L2ErrorPow2(function<double(RefPoint)> approximate, function<double(DomPoint)> exactSolution)
	{
		function<double(RefPoint)> errorFunction = [this, exactSolution, approximate](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return pow(exactSolution(domainPoint) - approximate(refElementPoint), 2);
		};

		CartesianShape<Dim>::ComputeIntegral(errorFunction, 0);
	}
};