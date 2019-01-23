#pragma once
#include "Element.h"
#include "Poisson_DG_ReferenceElement.h"

class Poisson_DG_ReferenceCube : public Poisson_DG_ReferenceElement
{
public:
	Poisson_DG_ReferenceCube(int nBasisFunctions) :
		Poisson_DG_ReferenceElement(nBasisFunctions)
	{}

	void ComputeVolumicTerm(BasisFunction* p_phi1, BasisFunction* p_phi2)
	{
		IBasisFunction3D* phi1 = static_cast<IBasisFunction3D*>(p_phi1);
		IBasisFunction3D* phi2 = static_cast<IBasisFunction3D*>(p_phi2);

		function<double(double, double, double)> functionToIntegrate = [phi1, phi2](double t, double u, double v) {
			return InnerProduct(phi1->Grad(t, u, v), phi2->Grad(t, u, v));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		double result = Utils::Integral(nQuadPoints, functionToIntegrate, -1,1, -1,1, -1,1);

		this->_volumicTerms(phi1->LocalNumber, phi2->LocalNumber) = result;
	}

	void ComputeMassTerm(BasisFunction* p_phi1, BasisFunction* p_phi2)
	{
		IBasisFunction3D* phi1 = static_cast<IBasisFunction3D*>(p_phi1);
		IBasisFunction3D* phi2 = static_cast<IBasisFunction3D*>(p_phi2);

		function<double(double, double, double)> functionToIntegrate = [phi1, phi2](double t, double u, double v) {
			return phi1->Eval(t, u, v)*phi2->Eval(t, u, v);
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		double result = Utils::Integral(nQuadPoints, functionToIntegrate, -1,1, -1,1, -1,1);

		this->_massTerms(phi1->LocalNumber, phi2->LocalNumber) = result;
	}

private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
	}
};