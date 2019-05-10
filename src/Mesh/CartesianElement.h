#pragma once
#include "Element.h"
#include "CartesianShape.h"

template <int Dim>
class CartesianElement : virtual public Element<Dim>, public CartesianShape<Dim>
{
public:

	CartesianElement(BigNumber number, Point origin, double width) : Element<Dim>(number), CartesianShape<Dim>(origin, width)
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
		Point origin = CartesianShape<Dim>::Origin;
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

	Point ConvertToDomain(Point referenceElementPoint)
	{
		return CartesianShape<Dim>::ConvertToDomain(referenceElementPoint);
	}

	Point ConvertToReference(Point domainPoint)
	{
		return CartesianShape<Dim>::ConvertToReference(domainPoint);
	}

	vector<Point> GetNodalPoints(FunctionalBasis<Dim>* basis)
	{
		return CartesianShape<Dim>::GetNodalPoints(basis);
	}

	double SourceTerm(BasisFunction<Dim>* phi, SourceFunction* f)
	{
		function<double(Point)> sourceTimesBasisFunction = [this, f, phi](Point refElementPoint) {
			Point domainPoint = this->ConvertToDomain(refElementPoint);
			return f->Eval(domainPoint) * phi->Eval(refElementPoint);
		};

		double h = this->Width;
		return pow(h / 2, Dim) *  Utils::Integral<Dim>(sourceTimesBasisFunction);
	}

	double L2ErrorPow2(function<double(Point)> approximate, function<double(Point)> exactSolution)
	{
		function<double(Point)> errorFunction = [this, exactSolution, approximate](Point refElementPoint) {
			Point domainPoint = this->ConvertToDomain(refElementPoint);
			return pow(exactSolution(domainPoint) - approximate(refElementPoint), 2);
		};

		double h = this->Width;
		return pow(h / 2, Dim) *  Utils::Integral<Dim>(errorFunction);
	}
};