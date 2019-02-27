#pragma once
#include "Face.h"
#include "../Utils/Utils.h"
#include "Cube.h"

class SquareFace : public Face<3>, public Poisson_DG_Face<3>
{
public:
	double Width = 0;

	SquareFace(BigNumber number, double width, Element<3>* element1, Element<3>* element2) : Face(number, element1, element2)
	{
		this->Width = width;
	}

	SquareFace(BigNumber number, double width, Element<3>* element1) : Face(number, element1)
	{
		this->Width = width;
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double CouplingTerm(Element<3>* element1, BasisFunction<3>* p_phi1, Element<3>* element2, BasisFunction<3>* p_phi2, DiffusionPartition diffusionPartition)
	{
		auto n1 = element1->OuterNormalVector(this);
		auto n2 = element2->OuterNormalVector(this);

		double k1 = element1->DiffusionCoefficient(diffusionPartition);
		double k2 = element2->DiffusionCoefficient(diffusionPartition);

		double weight1 = 1;
		double weight2 = 1;
		if (!this->IsDomainBoundary)
		{
			Element<3>* elementOnTheOtherSide1 = element1->ElementOnTheOtherSideOf(this);
			Element<3>* elementOnTheOtherSide2 = element2->ElementOnTheOtherSideOf(this);
			double l1 = k1;
			double l2 = elementOnTheOtherSide1->DiffusionCoefficient(diffusionPartition);
			weight1 = elementOnTheOtherSide1->DiffusionCoefficient(diffusionPartition) / (l1 + l2);
			weight2 = elementOnTheOtherSide2->DiffusionCoefficient(diffusionPartition) / (l1 + l2);
		}

		auto phi1 = element1->EvalPhiOnFace(this, p_phi1);
		auto gradPhi1 = element1->GradPhiOnFace(this, p_phi1);

		auto phi2 = element2->EvalPhiOnFace(this, p_phi2);
		auto gradPhi2 = element2->GradPhiOnFace(this, p_phi2);

		std::function<double(double, double)> functionToIntegrate = [n1, n2, phi1, phi2, gradPhi1, gradPhi2, weight1, weight2, k1, k2](double u, double v) {
			Point p(u, v);
			double meanGradPhi1_scal_jumpPhi2 = weight1 * k1 * InnerProduct(gradPhi1(p), n2) * phi2(p);
			double meanGradPhi2_scal_jumpPhi1 = weight2 * k2 * InnerProduct(gradPhi2(p), n1) * phi1(p);
			return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
		};

		int nQuadPoints = p_phi1->GetDegree() + p_phi2->GetDegree() + 1;
		double h = this->Width;
		return -h / 2 * Utils::Integral(nQuadPoints, functionToIntegrate, -1,1, -1,1);
	}

	double PenalizationTerm(Poisson_DG_Element<3>* element1, BasisFunction<3>* p_phi1, Poisson_DG_Element<3>* element2, BasisFunction<3>* p_phi2, double penalizationCoefficient, DiffusionPartition diffusionPartition)
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

		double diffusionDependantCoefficient = element1->DiffusionCoefficient(diffusionPartition);
		if (!this->IsDomainBoundary)
		{
			double k1 = this->Element1->DiffusionCoefficient(diffusionPartition);
			double k2 = this->Element2->DiffusionCoefficient(diffusionPartition);
			diffusionDependantCoefficient = 2 * k1*k2 / (k1 + k2);
		}

		return diffusionDependantCoefficient * penalizationCoefficient * integralJump1ScalarJump2;
	}

private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
	}
};