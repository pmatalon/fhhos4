#pragma once
#include "Face.h"
#include "../Utils/Utils.h"
#include "../DG/Poisson_DG_Face.h"

class IntervalFace : public Face, public Poisson_DG_Face<2>
{
public:
	double Length;

	IntervalFace(BigNumber number, double length, Element* element1, Element* element2) : Face(number, element1, element2)
	{	
		this->Length = length;
	}

	IntervalFace(BigNumber number, double length, Element* element1) : Face(number, element1)
	{	
		this->Length = length;
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double CouplingTerm(Poisson_DG_Element<2>* element1, BasisFunction<2>* p_phi1, Poisson_DG_Element<2>* element2, BasisFunction<2>* p_phi2)
	{
		auto n1 = element1->OuterNormalVector(this);
		auto n2 = element2->OuterNormalVector(this);

		double meanFactor = this->IsDomainBoundary ? 1 : 0.5;

		auto phi1 = element1->EvalPhiOnFace(this, p_phi1);
		auto gradPhi1 = element1->GradPhiOnFace(this, p_phi1);

		auto phi2 = element2->EvalPhiOnFace(this, p_phi2);
		auto gradPhi2 = element2->GradPhiOnFace(this, p_phi2);

		std::function<double(double)> functionToIntegrate = [n1, n2, phi1, phi2, gradPhi1, gradPhi2](double u) {
			Point p(u);
			double meanGradPhi1_scal_jumpPhi2 = InnerProduct(gradPhi1(p), n2) * phi2(p);
			double meanGradPhi2_scal_jumpPhi1 = InnerProduct(gradPhi2(p), n1) * phi1(p);
			return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
		};

		int nQuadPoints = p_phi1->GetDegree() + p_phi2->GetDegree() + 1;
		return -meanFactor * Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1);
	}

	double PenalizationTerm(Poisson_DG_Element<2>* element1, BasisFunction<2>* p_phi1, Poisson_DG_Element<2>* element2, BasisFunction<2>* p_phi2, double penalizationCoefficient)
	{
		auto n1 = element1->OuterNormalVector(this);
		auto n2 = element2->OuterNormalVector(this);

		assert(InnerProduct(n1, n2) == 1 || InnerProduct(n1, n2) == -1);

		auto phi1 = element1->EvalPhiOnFace(this, p_phi1);
		auto phi2 = element2->EvalPhiOnFace(this, p_phi2);

		std::function<double(double)> functionToIntegrate = [phi1, phi2, n1, n2](double s) {
			Point p(s);
			return InnerProduct(n1, n2) * phi1(p) * phi2(p);
		};

		int nQuadPoints = p_phi1->GetDegree() + p_phi2->GetDegree() + 2;
		double h = this->Length;
		double integralJump1ScalarJump2 = h / 2 * Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1);
		return penalizationCoefficient * integralJump1ScalarJump2;
	}

private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0] + vector1[1] * vector2[1];
	}
};