#pragma once
#include "../Mesh/Element.h"
#include "../FunctionalBasis/BasisFunction.h"

template <int Dim>
class Poisson_DG_Face : virtual public Face<Dim>
{
public:
	Poisson_DG_Face(BigNumber number, Element<Dim>* element1, Element<Dim>* element2) : Face<Dim>(number, element1, element2) {}

	virtual double CouplingTerm(Element<Dim>* element1, BasisFunction<Dim>* p_phi1, Element<Dim>* element2, BasisFunction<Dim>* p_phi2)
	{
		auto n1 = element1->OuterNormalVector(this);
		auto n2 = element2->OuterNormalVector(this);

		double k1 = element1->Kappa;
		double k2 = element2->Kappa;

		double weight1 = 1;
		double weight2 = 1;
		if (!this->IsDomainBoundary)
		{
			Element<Dim>* elementOnTheOtherSide1 = element1->ElementOnTheOtherSideOf(this);
			Element<Dim>* elementOnTheOtherSide2 = element2->ElementOnTheOtherSideOf(this);
			double l1 = k1;
			double l2 = elementOnTheOtherSide1->Kappa;
			weight1 = elementOnTheOtherSide1->Kappa / (l1 + l2);
			weight2 = elementOnTheOtherSide2->Kappa / (l1 + l2);
		}

		auto phi1 = element1->EvalPhiOnFace(this, p_phi1);
		auto gradPhi1 = element1->GradPhiOnFace(this, p_phi1);

		auto phi2 = element2->EvalPhiOnFace(this, p_phi2);
		auto gradPhi2 = element2->GradPhiOnFace(this, p_phi2);

		RefFunction functionToIntegrate = [n1, n2, phi1, phi2, gradPhi1, gradPhi2, weight1, weight2, k1, k2](RefPoint p) {
			double meanGradPhi1_scal_jumpPhi2 = weight1 * k1 * gradPhi1(p).dot(n2) * phi2(p);
			double meanGradPhi2_scal_jumpPhi1 = weight2 * k2 * gradPhi2(p).dot(n1) * phi1(p);
			return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
		};

		int polynomialDegree = max(p_phi1->GetDegree() + p_phi2->GetDegree() - 1, 0);
		double integralJump1ScalarJump2 = this->Integral(functionToIntegrate, polynomialDegree);
		return -integralJump1ScalarJump2;
	}

	virtual double PenalizationTerm(Element<Dim>* element1, BasisFunction<Dim>* p_phi1, Element<Dim>* element2, BasisFunction<Dim>* p_phi2, double penalizationCoefficient)
	{
		auto n1 = element1->OuterNormalVector(this);
		auto n2 = element2->OuterNormalVector(this);

		auto phi1 = element1->EvalPhiOnFace(this, p_phi1);
		auto phi2 = element2->EvalPhiOnFace(this, p_phi2);

		RefFunction functionToIntegrate = [phi1, phi2, n1, n2](RefPoint p) {
			return n1.dot(n2) * phi1(p) * phi2(p);
		};

		int polynomialDegree = p_phi1->GetDegree() + p_phi2->GetDegree();
		double integralJump1ScalarJump2 = this->Integral(functionToIntegrate, polynomialDegree);

		double diffusionDependantCoefficient = element1->Kappa;
		if (!this->IsDomainBoundary)
		{
			double k1 = this->Element1->Kappa;
			double k2 = this->Element2->Kappa;
			diffusionDependantCoefficient = 2 * k1*k2 / (k1 + k2);
		}

		return diffusionDependantCoefficient * penalizationCoefficient * integralJump1ScalarJump2;
	}
};