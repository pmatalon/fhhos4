#pragma once
#include "Element.h"
#include "Poisson_DG_ReferenceElement.h"

class Poisson_DG_ReferenceInterval : public Poisson_DG_ReferenceElement
{
public:
	// [-1, 1]
	Poisson_DG_ReferenceInterval(int nBasisFunctions) : 
		Poisson_DG_ReferenceElement(nBasisFunctions)
	{}

	void ComputeVolumicTerm(BasisFunction* p_phi1, BasisFunction* p_phi2)
	{
		IBasisFunction1D* phi1 = dynamic_cast<IBasisFunction1D*>(p_phi1);
		IBasisFunction1D* phi2 = dynamic_cast<IBasisFunction1D*>(p_phi2);

		DefInterval refInterval = phi1->DefinitionInterval();

		function<double(double)> functionToIntegrate = [phi1, phi2](double t) {
			return InnerProduct(phi1->Grad(t), phi2->Grad(t));
		};
		// Sans doute probleme ici pour Bernstein
		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		double factor = refInterval.Length / 2;
		double result = factor * Utils::Integral(nQuadPoints, functionToIntegrate, refInterval);

		this->_volumicTerms(phi1->LocalNumber, phi2->LocalNumber) = result;
	}

private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0];
	}
};