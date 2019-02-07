#pragma once
#include <math.h>
#include "Poisson_DGTerms.h"
#include "BasisFunction.h"
#include "Utils.h"
#include "Element.h"
#include "Face.h"
#include "Poisson_DG_ReferenceSquare.h"
#include "Square.h"
using namespace std;

class Poisson2D_DGTerms_LocalBasis : public Poisson_DGTerms
{
public:
	Poisson2D_DGTerms_LocalBasis(SourceFunction* sourceFunction, FunctionalBasis2D* basis)
		: Poisson_DGTerms(sourceFunction)
	{
		Poisson_DG_ReferenceElement* refSquare = new Poisson_DG_ReferenceSquare(basis->NumberOfLocalFunctionsInElement(NULL));
		this->ComputeReferenceTerms(basis, refSquare);
		this->ReferenceElements.insert(std::make_pair(StandardElementCode::Square, refSquare));
	}

	bool IsGlobalBasis() { return false; }

	/*double VolumicTerm(Element* element, IBasisFunction2D* phi1, IBasisFunction2D* phi2)
	{
		Square* square = static_cast<Square*>(element);
		return this->VolumicTerm(square, phi1, phi2);
	}*/

	/*double VolumicTerm(Square* element, IBasisFunction2D* phi1, IBasisFunction2D* phi2)
	{
		function<double(double, double)> functionToIntegrate = [phi1, phi2](double t, double u) {
			return InnerProduct(phi1->GradPhiOnFace(t, u), phi2->GradPhiOnFace(t, u));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		return Utils::Integral(nQuadPoints, functionToIntegrate, phi1->DefinitionInterval(), phi1->DefinitionInterval());
	}*/

	/*double CouplingTerm(Face* face, Element* element1, IBasisFunction2D* phi1, Element* element2, IBasisFunction2D* phi2)
	{
		assert(face->IsBetween(element1, element2));

		Square* square1 = static_cast<Square*>(element1);
		Square* square2 = static_cast<Square*>(element2);
		return this->CouplingTerm(face, square1, phi1, square2, phi2);
	}

	double CouplingTerm(Face* face, Square* element1, IBasisFunction2D* phi1, Square* element2, IBasisFunction2D* phi2)
	{
		auto n1 = element1->OuterNormalVector(face);
		auto n2 = element2->OuterNormalVector(face);

		double meanFactor = face->IsDomainBoundary ? 1 : 0.5;

		Face2D* face2D = (Face2D*)face;

		std::function<double(double)> functionToIntegrate;
		if (face2D->IsVertical())
		{
			double t1 = face == element1->EastFace ? 1 : -1;
			double t2 = face == element2->EastFace ? 1 : -1;

			functionToIntegrate = [n1, n2, phi1, phi2, t1, t2](double u) {
				double meanGradPhi1_scal_jumpPhi2 = InnerProduct(phi1->GradPhiOnFace(t1, u), n2) * phi2->EvalPhiOnFace(t2, u);
				double meanGradPhi2_scal_jumpPhi1 = InnerProduct(phi2->GradPhiOnFace(t2, u), n1) * phi1->EvalPhiOnFace(t1, u);
				return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
			};
		}
		else if (face2D->IsHorizontal())
		{
			double u1 = face == element1->NorthFace ? 1 : -1;
			double u2 = face == element2->NorthFace ? 1 : -1;

			functionToIntegrate = [n1, n2, phi1, phi2, u1, u2](double t) {
				double meanGradPhi1_scal_jumpPhi2 = InnerProduct(phi1->GradPhiOnFace(t, u1), n2) * phi2->EvalPhiOnFace(t, u2);
				double meanGradPhi2_scal_jumpPhi1 = InnerProduct(phi2->GradPhiOnFace(t, u2), n1) * phi1->EvalPhiOnFace(t, u1);
				return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
			};
		}
		else
			assert(false);

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 1;
		return -meanFactor * Utils::Integral(nQuadPoints, functionToIntegrate, -1,1);
	}*/

	/*double PenalizationTerm(Face* face, Element* element1, IBasisFunction2D* phi1, Element* element2, IBasisFunction2D* phi2, double penalizationCoefficient) override
	{
		assert(face->IsBetween(element1, element2));

		Square* square1 = static_cast<Square*>(element1);
		Square* square2 = static_cast<Square*>(element2);
		return this->PenalizationTerm(face, square1, phi1, square2, phi2, penalizationCoefficient);
	}

	double PenalizationTerm(Face* face, Square* element1, IBasisFunction2D* phi1, Square* element2, IBasisFunction2D* phi2, double penalizationCoefficient)
	{
		double h = element1->Width;
		auto n1 = element1->OuterNormalVector(face);
		auto n2 = element2->OuterNormalVector(face);

		Face2D* face2D = (Face2D*)face;

		function<double(double)> functionToIntegrate;

		if (face2D->IsVertical())
		{
			double t1 = face == element1->EastFace ? 1 : -1;
			double t2 = face == element2->EastFace ? 1 : -1;

			functionToIntegrate = [phi1, phi2, n1, n2, t1, t2](double u) {
				return InnerProduct(n1, n2) * phi1->Eval(t1, u) * phi2->Eval(t2, u);
			};
		}
		else if (face2D->IsHorizontal())
		{
			double u1 = face == element1->NorthFace ? 1 : -1;
			double u2 = face == element2->NorthFace ? 1 : -1;

			functionToIntegrate = [phi1, phi2, n1, n2, u1, u2](double t) {
				return InnerProduct(n1, n2) * phi1->Eval(t, u1) * phi2->Eval(t, u2);
			};
		}
		else
			assert(false);

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		double jacobian = h / 2;
		double integralJump1ScalarJump2 = jacobian * Utils::Integral(nQuadPoints, functionToIntegrate, -1,1);
		return penalizationCoefficient * integralJump1ScalarJump2;
	}*/

	/*double RightHandSide(Element* element, IBasisFunction2D* phi)
	{
		Square* square = static_cast<Square*>(element);
		return this->RightHandSide(square, phi);
	}

	double RightHandSide(Square* element, IBasisFunction2D* phi)
	{
		double x1 = element->X;
		double x2 = element->X + element->Width;
		double y1 = element->Y;
		double y2 = element->Y + element->Width;

		function<double(double, double)> sourceTimesBasisFunction = [this, phi, x1, x2, y1, y2](double t, double u) {
			return this->_sourceFunction((x2 - x1) / 2 * t + (x2 + x1) / 2, (y2 - y1) / 2 * u + (y2 + y1) / 2) * phi->Eval(t, u);
		};

		double jacobian = (x2 - x1) * (y2 - y1) / 4;
		return jacobian * Utils::Integral(sourceTimesBasisFunction, -1,1, -1,1);
	}*/

/*private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0] + vector1[1] * vector2[1];
	}*/
};