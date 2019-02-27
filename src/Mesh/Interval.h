#pragma once
#include "Element.h"
#include "../DG/Poisson_DG_Element.h"
#include "../DG/Poisson_DG_ReferenceElement.h"

class Interval : public Element<1>, public Poisson_DG_Element<1>
{
public:
	double A;
	double B;

	Face<1>* Left;
	Face<1>* Right;

	Interval(BigNumber number, double a, double b, Face<1>* left, Face<1>* right) : Element(number)
	{
		this->A = a;
		this->B = b;
		this->Faces.push_back(left);
		this->Faces.push_back(right);
		this->Left = left;
		this->Right = right;
	}

	StandardElementCode StdElementCode()
	{
		return StandardElementCode::Interval;
	}

	double* OuterNormalVector(Face<1>* interface)
	{
		if (interface == this->Left)
			return new double[1]{ -1 };
		else if(interface == this->Right)
			return new double[1]{ 1 };
		return NULL;
	}

	double DiffusionCoefficient(DiffusionPartition diffusionPartition)
	{
		return diffusionPartition.Coefficient(Point((this->A + this->B) / 2));
	}

	double Integral(function<double(Point)> func)
	{
		function<double(double)> funcToIntegrate = [func](double x) {
			return func(x);
		};

		return Utils::Integral(funcToIntegrate, this->A, this->B);
	}

	double L2ErrorPow2(function<double(Point)> approximate, function<double(Point)> exactSolution)
	{
		double a = this->A;
		double b = this->B;

		function<double(double)> errorFunction = [exactSolution, approximate, a, b](double t) {
			return pow(exactSolution((b - a) / 2 * t + (b + a) / 2) - approximate(t), 2);
		};

		return (b - a) / 2 * Utils::Integral(errorFunction, -1, 1);
	}
	
	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double VolumicTerm(BasisFunction<1>* phi1, BasisFunction<1>* phi2, Poisson_DG_ReferenceElement<1>* referenceElement, DiffusionPartition diffusionPartition)
	{
		double h = this->B - this->A;
		double kappa = this->DiffusionCoefficient(diffusionPartition);
		return 2 / h * kappa * referenceElement->VolumicTerm(phi1, phi2);
	}

	double MassTerm(BasisFunction<1>* phi1, BasisFunction<1>* phi2, Poisson_DG_ReferenceElement<1>* referenceElement)
	{
		double h = this->B - this->A;
		return h / 2 * referenceElement->MassTerm(phi1, phi2);
	}
	
	double SourceTerm(BasisFunction<1>* phi, SourceFunction* f)
	{
		double a = this->A;
		double b = this->B;

		function<double(double)> sourceTimesBasisFunction = [f, phi, a, b](double t) {
			return f->Eval(Point((b - a) / 2 * t + (a + b) / 2)) * phi->Eval(Point(t));
		};

		return  (b - a) / 2 * Utils::Integral(sourceTimesBasisFunction, -1, 1);
	}

	function<double(Point)> EvalPhiOnFace(Face<1>* face, BasisFunction<1>* p_phi)
	{
		IBasisFunction1D* phi = dynamic_cast<IBasisFunction1D*>(p_phi);

		double tFixed = face == this->Left ? -1 : 1;
		function<double(Point)> evalOnFace = [phi, tFixed](Point point0D) {
			return phi->Eval(tFixed);
		};
		return evalOnFace;
	}


	function<double*(Point)> GradPhiOnFace(Face<1>* face, BasisFunction<1>* p_phi)
	{
		IBasisFunction1D* phi = dynamic_cast<IBasisFunction1D*>(p_phi);

		double tFixed = face == this->Left ? -1 : 1;
		function<double*(Point)> gradOnFace = [phi, tFixed](Point point0D) {
			return phi->Grad(tFixed);
		};
		return gradOnFace;
	}
};