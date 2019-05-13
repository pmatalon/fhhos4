#pragma once
#include "CartesianElement.h"
#include "../DG/Poisson_DG_Element.h"

class Interval : public CartesianElement<1>, public Poisson_DG_Element<1>
{
public:
	Face<1>* Left;
	Face<1>* Right;

	Interval(BigNumber number, double a, double b, Face<1>* left, Face<1>* right) : Element(number), CartesianElement(number, DomPoint(a), b-a), Poisson_DG_Element(number)
	{
		this->AddFace(left);
		this->AddFace(right);
		this->Left = left;
		this->Right = right;
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	vector<double> OuterNormalVector(Face<1>* interface)
	{
		if (interface == this->Left)
			return vector<double>{ -1 };
		else if(interface == this->Right)
			return vector<double>{ 1 };
		assert(false);
	}

	double IntegralGlobalFunction(function<double(DomPoint)> func)
	{
		function<double(double)> funcToIntegrate = [func](double x) {
			return func(x);
		};

		return Utils::Integral(funcToIntegrate, this->Origin.X, this->Origin.X + this->Width);
	}

	
	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double VolumicTerm(BasisFunction<1>* phi1, BasisFunction<1>* phi2, DiffusionPartition diffusionPartition)
	{
		double kappa = CartesianElement::DiffusionCoefficient(diffusionPartition);
		return kappa * CartesianElement::IntegralGradGrad(phi1, phi2);
	}

	double MassTerm(BasisFunction<1>* phi1, BasisFunction<1>* phi2)
	{
		return CartesianElement::MassTerm(phi1, phi2);
	}
	
	double SourceTerm(BasisFunction<1>* phi, SourceFunction* f)
	{
		return CartesianElement::SourceTerm(phi, f);
	}
};