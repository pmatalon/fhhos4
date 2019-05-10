#pragma once
#include <Eigen/Dense>
#include "../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
class ReferenceCartesianShape
{
private:
	static ReferenceCartesianShape* _instance;

	Eigen::MatrixXd _stiffnessMatrix;
	Eigen::MatrixXd _massMatrix;
public:
	static ReferenceCartesianShape* Get()
	{
		return _instance;
	}

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->_massMatrix(phi1->LocalNumber, phi2->LocalNumber);
	}
	double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->_stiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber);
	}


private:

	ReferenceCartesianShape()
	{
	}

	double Integral(BasisFunction<Dim>* phi)
	{
		return Utils::Integral(phi);
	}

	double ComputeMassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		function<double(Point)> functionToIntegrate = [phi1, phi2](RefPoint p) {
			return phi1->Eval(p)*phi2->Eval(p);
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		return Utils::Integral<Dim>(nQuadPoints, functionToIntegrate);
	}

	double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi1->GetDegree() == 0)
			return 0;

		function<double(Point)> functionToIntegrate = [phi1, phi2](RefPoint p) {
			return Element<Dim>::InnerProduct(phi1->Grad(p), phi2->Grad(p));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		return Utils::Integral<Dim>(nQuadPoints, functionToIntegrate);
	}
};