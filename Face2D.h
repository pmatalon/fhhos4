#pragma once
#include "Face.h"
#include "Utils.h"
#include "Poisson_DG_Face.h"

class Face2D : public Face, public Poisson_DG_Face
{
public:
	double Length;
	double X1;
	double Y1;
	double X2;
	double Y2;

	Face2D(BigNumber number, double length, Element* element1, Element* element2) : Face(number, element1, element2)
	{	
		this->Length = length;
	}

	Face2D(BigNumber number, double length, Element* element1) : Face(number, element1)
	{	
		this->Length = length;
	}

	bool IsVertical() {	return this->X1 == this->X2; }
	bool IsHorizontal()	{ return this->Y1 == this->Y2; }

	string ToString() override
	{
		string s = "Interface " + std::to_string(this->Number);
		if (this->IsVertical())
			s += " (vertical)";
		else if (this->IsHorizontal())
			s += " (horizontal)";
		if (this->IsDomainBoundary)
			s += " on element " + std::to_string(this->Element1->Number) + " (boundary)";
		else
			s += " between element " + std::to_string(this->Element1->Number) + " and element " + std::to_string(this->Element2->Number);
		return s;
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double CouplingTerm(Poisson_DG_Element* element1, BasisFunction* p_phi1, Poisson_DG_Element* element2, BasisFunction* p_phi2)
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

	double PenalizationTerm(Poisson_DG_Element* element1, BasisFunction* p_phi1, Poisson_DG_Element* element2, BasisFunction* p_phi2, double penalizationCoefficient)
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