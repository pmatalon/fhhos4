#pragma once
#include "Element.h"

template <short Dim>
class CartesianElement : public Element<Dim>
{
public:
	double Width;

	CartesianElement(int number, double width) : Element<Dim>(number)
	{
		this->Width = width;
	}

	double GetDiameter() override
	{
		return this->Width;
	}

	/*double Integral(function<double(Point)> func) override
	{

	}*/

	double Integral(BasisFunction<Dim>* phi)
	{
		double h = this->Width;
		return pow(h/2, Dim) * Utils::Integral(phi);
	}

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) override
	{
		double h = this->Width;

		function<double(Point)> functionToIntegrate = [phi1, phi2](Point p) {
			return phi1->Eval(p)*phi2->Eval(p);
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		return pow(h / 2, Dim) * Utils::Integral<Dim>(nQuadPoints, functionToIntegrate);
	}

	double IntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi1->GetDegree() == 0)
			return 0;

		double h = this->Width;

		function<double(Point)> functionToIntegrate = [phi1, phi2](Point p) {
			return Element<Dim>::InnerProduct(phi1->Grad(p), phi2->Grad(p));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		return pow(h / 2, Dim-2) * Utils::Integral<Dim>(nQuadPoints, functionToIntegrate);
	}

	virtual Point ConvertToDomain(Point referenceElementPoint) = 0;

	virtual double SourceTerm(BasisFunction<Dim>* phi, SourceFunction* f)
	{
		function<double(Point)> sourceTimesBasisFunction = [this, f, phi](Point refElementPoint) {
			Point domainPoint = this->ConvertToDomain(refElementPoint);
			return f->Eval(domainPoint) * phi->Eval(refElementPoint);
		};

		double h = this->Width;
		return pow(h / 2, Dim) *  Utils::Integral<Dim>(sourceTimesBasisFunction);
	}

	virtual double L2ErrorPow2(function<double(Point)> approximate, function<double(Point)> exactSolution)
	{
		function<double(Point)> errorFunction = [this, exactSolution, approximate](Point refElementPoint) {
			Point domainPoint = this->ConvertToDomain(refElementPoint);
			return pow(exactSolution(domainPoint) - approximate(refElementPoint), 2);
		};

		double h = this->Width;
		return pow(h / 2, Dim) *  Utils::Integral<Dim>(errorFunction);
	}
};