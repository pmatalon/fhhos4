#pragma once
#include "Face.h"
#include "Utils.h"
#include "Cube.h"

class Face3D : public Face, public Poisson_DG_Face<3>
{
public:
	double Width = 0;

	Face3D(BigNumber number, double width, Element* element1, Element* element2) : Face(number, element1, element2)
	{
		this->Width = width;
	}

	Face3D(BigNumber number, double width, Element* element1) : Face(number, element1)
	{
		this->Width = width;
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double CouplingTerm(Poisson_DG_Element<3>* element1, BasisFunction<3>* p_phi1, Poisson_DG_Element<3>* element2, BasisFunction<3>* p_phi2)
	{
		auto n1 = element1->OuterNormalVector(this);
		auto n2 = element2->OuterNormalVector(this);

		double meanFactor = this->IsDomainBoundary ? 1 : 0.5;

		auto phi1 = element1->EvalPhiOnFace(this, p_phi1);
		auto gradPhi1 = element1->GradPhiOnFace(this, p_phi1);

		auto phi2 = element2->EvalPhiOnFace(this, p_phi2);
		auto gradPhi2 = element2->GradPhiOnFace(this, p_phi2);

		std::function<double(double, double)> functionToIntegrate = [n1, n2, phi1, phi2, gradPhi1, gradPhi2](double u, double v) {
			Point p(u, v);
			double meanGradPhi1_scal_jumpPhi2 = InnerProduct(gradPhi1(p), n2) * phi2(p);
			double meanGradPhi2_scal_jumpPhi1 = InnerProduct(gradPhi2(p), n1) * phi1(p);
			return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
		};

		int nQuadPoints = p_phi1->GetDegree() + p_phi2->GetDegree() + 1;
		double h = this->Width;
		return -meanFactor * h / 2 * Utils::Integral(nQuadPoints, functionToIntegrate, -1,1, -1,1);
	}

	double PenalizationTerm(Poisson_DG_Element<3>* element1, BasisFunction<3>* p_phi1, Poisson_DG_Element<3>* element2, BasisFunction<3>* p_phi2, double penalizationCoefficient)
	{
		auto n1 = element1->OuterNormalVector(this);
		auto n2 = element2->OuterNormalVector(this);

		auto phi1 = element1->EvalPhiOnFace(this, p_phi1);
		auto phi2 = element2->EvalPhiOnFace(this, p_phi2);

		std::function<double(double, double)> functionToIntegrate = [phi1, phi2, n1, n2](double s, double t) {
			Point p(s, t);
			return InnerProduct(n1, n2) * phi1(p) * phi2(p);
		};

		int nQuadPoints = p_phi1->GetDegree() + p_phi2->GetDegree() + 2;
		double h = this->Width;
		double integralJump1ScalarJump2 = pow(h, 2) / 4 * Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1, -1, 1);
		return penalizationCoefficient * integralJump1ScalarJump2;
	}

private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
	}
};