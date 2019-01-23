#pragma once
#include "Element.h"
#include "Poisson_DG_ReferenceElement.h"

class Poisson_DG_ReferenceSquare : public Poisson_DG_ReferenceElement
{
public:
	Poisson_DG_ReferenceSquare(int nBasisFunctions) :
		Poisson_DG_ReferenceElement(nBasisFunctions)
	{}

	void ComputeVolumicTerm(BasisFunction* p_phi1, BasisFunction* p_phi2)
	{
		IBasisFunction2D* phi1 = dynamic_cast<IBasisFunction2D*>(p_phi1);
		IBasisFunction2D* phi2 = dynamic_cast<IBasisFunction2D*>(p_phi2);

		function<double(double, double)> functionToIntegrate = [phi1, phi2](double t, double u) {
			return InnerProduct(phi1->Grad(t, u), phi2->Grad(t, u));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		double result = Utils::Integral(nQuadPoints, functionToIntegrate, -1,1, -1,1);

		this->_volumicTerms(phi1->LocalNumber, phi2->LocalNumber) = result;
	}

	void ComputeMassTerm(BasisFunction* p_phi1, BasisFunction* p_phi2)
	{
		IBasisFunction2D* phi1 = dynamic_cast<IBasisFunction2D*>(p_phi1);
		IBasisFunction2D* phi2 = dynamic_cast<IBasisFunction2D*>(p_phi2);

		function<double(double, double)> functionToIntegrate = [phi1, phi2](double t, double u) {
			return phi1->Eval(t, u)*phi2->Eval(t, u);
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		double result = Utils::Integral(nQuadPoints, functionToIntegrate, -1,1, -1,1);

		this->_massTerms(phi1->LocalNumber, phi2->LocalNumber) = result;
	}

private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0] + vector1[1] * vector2[1];
	}
};