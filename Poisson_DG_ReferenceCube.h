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
		IBasisFunction3D* phi1 = dynamic_cast<IBasisFunction3D*>(p_phi1);
		IBasisFunction3D* phi2 = dynamic_cast<IBasisFunction3D*>(p_phi2);

		function<double(double, double, double)> functionToIntegrate = [phi1, phi2](double t, double u, double v) {
			return InnerProduct(phi1->Grad(t, u, v), phi2->Grad(t, u, v));
		};

		DefInterval refInterval = phi1->DefinitionInterval();
		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		double result = Utils::Integral(nQuadPoints, functionToIntegrate, refInterval, refInterval, refInterval);

		this->_volumicTerms(phi1->LocalNumber, phi2->LocalNumber) = result;
	}

private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
	}
};