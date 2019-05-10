#pragma once
#include "Face.h"
#include "../Utils/Utils.h"
#include "../DG/Poisson_DG_Face.h"
#include "../HHO/Poisson_HHO_Face.h"

class IntervalFace : virtual public Face<2>, public CartesianShape<2,1>, public Poisson_DG_Face<2>, public Poisson_HHO_Face<2>
{
public:

	IntervalFace(BigNumber number, Point origin, double length, Element<2>* element1, Element<2>* element2, CartesianShapeOrientation orientation) : 
		Face(number, element1, element2), 
		CartesianShape(origin, length, orientation),
		Poisson_DG_Face(number, element1, element2),
		Poisson_HHO_Face(number, element1, element2)
	{
	}

	IntervalFace(BigNumber number, Point origin, double length, Element<2>* element1, CartesianShapeOrientation orientation) :
		Face(number, element1), 
		CartesianShape(origin, length, orientation),
		Poisson_DG_Face(number, element1, NULL),
		Poisson_HHO_Face(number, element1, NULL)
	{
	}

	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	double GetDiameter()
	{
		return CartesianShape::Width;
	}

	double Measure()
	{
		return CartesianShape::Measure();
	}

	double MassTerm(BasisFunction<1>* phi1, BasisFunction<1>* phi2) override
	{
		return CartesianShape::MassTerm(phi1, phi2);
	}

	double MassTerm(BasisFunction<1>* facePhi, Element<2>* element, BasisFunction<2>* reconstructPhi) override
	{
		double h = this->Width;
		auto reconstructPhiOnFace = element->EvalPhiOnFace(this, reconstructPhi);

		function<double(double)> functionToIntegrate = [facePhi, reconstructPhiOnFace](double t) {
			return facePhi->Eval(t) * reconstructPhiOnFace(t);
		};

		int nQuadPoints = facePhi->GetDegree() + reconstructPhi->GetDegree() + 2;
		return h / 2 * Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1);
	}

	Point ConvertToReference(Point domainPoint)
	{
		return CartesianShape<2,1>::ConvertToReference(domainPoint);
	}

	Point ConvertToDomain(Point referenceElementPoint)
	{
		return CartesianShape<2,1>::ConvertToDomain(referenceElementPoint);
	}

	vector<Point> GetNodalPoints(FunctionalBasis<1>* basis)
	{
		return CartesianShape<2,1>::GetNodalPoints(basis);
	}

	//---------------------------------------------------------------//
	//                 Poisson_DG_Face implementation                //
	//---------------------------------------------------------------//

	double CouplingTerm(Element<2>* element1, BasisFunction<2>* p_phi1, Element<2>* element2, BasisFunction<2>* p_phi2, DiffusionPartition diffusionPartition)
	{
		auto n1 = element1->OuterNormalVector(this);
		auto n2 = element2->OuterNormalVector(this);

		double k1 = element1->DiffusionCoefficient(diffusionPartition);
		double k2 = element2->DiffusionCoefficient(diffusionPartition);

		double weight1 = 1;
		double weight2 = 1;
		if (!this->IsDomainBoundary)
		{
			Element<2>* elementOnTheOtherSide1 = element1->ElementOnTheOtherSideOf(this);
			Element<2>* elementOnTheOtherSide2 = element2->ElementOnTheOtherSideOf(this);
			double l1 = k1;
			double l2 = elementOnTheOtherSide1->DiffusionCoefficient(diffusionPartition);
			weight1 = elementOnTheOtherSide1->DiffusionCoefficient(diffusionPartition) / (l1 + l2);
			weight2 = elementOnTheOtherSide2->DiffusionCoefficient(diffusionPartition) / (l1 + l2);
		}

		auto phi1 = element1->EvalPhiOnFace(this, p_phi1);
		auto gradPhi1 = element1->GradPhiOnFace(this, p_phi1);

		auto phi2 = element2->EvalPhiOnFace(this, p_phi2);
		auto gradPhi2 = element2->GradPhiOnFace(this, p_phi2);

		std::function<double(double)> functionToIntegrate = [n1, n2, phi1, phi2, gradPhi1, gradPhi2, weight1, weight2, k1, k2](double u) {
			Point p(u);
			double meanGradPhi1_scal_jumpPhi2 = weight1 * k1 * InnerProduct(gradPhi1(p), n2) * phi2(p);
			double meanGradPhi2_scal_jumpPhi1 = weight2 * k2 * InnerProduct(gradPhi2(p), n1) * phi1(p);
			return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
		};

		int nQuadPoints = p_phi1->GetDegree() + p_phi2->GetDegree() + 1;
		return -Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1);
	}

	double PenalizationTerm(Element<2>* element1, BasisFunction<2>* p_phi1, Element<2>* element2, BasisFunction<2>* p_phi2, double penalizationCoefficient, DiffusionPartition diffusionPartition)
	{
		auto n1 = element1->OuterNormalVector(this);
		auto n2 = element2->OuterNormalVector(this);

		auto phi1 = element1->EvalPhiOnFace(this, p_phi1);
		auto phi2 = element2->EvalPhiOnFace(this, p_phi2);

		std::function<double(double)> functionToIntegrate = [phi1, phi2, n1, n2](double s) {
			Point p(s);
			return InnerProduct(n1, n2) * phi1(p) * phi2(p);
		};

		int nQuadPoints = p_phi1->GetDegree() + p_phi2->GetDegree() + 2;
		double h = this->Width;
		double integralJump1ScalarJump2 = h / 2 * Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1);

		double diffusionDependantCoefficient = element1->DiffusionCoefficient(diffusionPartition);
		if (!this->IsDomainBoundary)
		{
			double k1 = this->Element1->DiffusionCoefficient(diffusionPartition);
			double k2 = this->Element2->DiffusionCoefficient(diffusionPartition);
			diffusionDependantCoefficient = 2 * k1*k2 / (k1 + k2);
		}

		return diffusionDependantCoefficient * penalizationCoefficient * integralJump1ScalarJump2;
	}

	friend ostream& operator<<(ostream& os, const IntervalFace& face)
	{
		os << face.Number << "\t(" << face.CartesianShape::Origin.X << ", " << face.CartesianShape::Origin.Y << ")\t" << (face.CartesianShape::Orientation == CartesianShapeOrientation::Horizontal ? "horizontal" : "vertical");
		return os;
	}

private:
	static double InnerProduct(vector<double> vector1, vector<double> vector2)
	{
		return vector1[0] * vector2[0] + vector1[1] * vector2[1];
	}
};