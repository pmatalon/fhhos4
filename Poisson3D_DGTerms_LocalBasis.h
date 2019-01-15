#pragma once
#include <functional>
#include <math.h>
#include "IPoisson_DGTerms.h"
#include "IBasisFunction.h"
#include "Utils.h"
#include "Element.h"
#include "Face.h"
#include "Cube.h"
#include "Poisson_DG_ReferenceCube.h"
using namespace std;

class Poisson3D_DGTerms_LocalBasis : public IPoisson_DGTerms<IBasisFunction3D>
{
protected:
	function<double(double, double, double)> _sourceFunction;

public:
	Poisson3D_DGTerms_LocalBasis(function<double(double, double, double)> sourceFunction, FunctionalBasis3D* basis)
	{
		this->_sourceFunction = sourceFunction;
		Poisson_DG_ReferenceElement* refCube = new Poisson_DG_ReferenceCube(basis->NumberOfLocalFunctionsInElement(NULL));
		this->ComputeVolumicTerms(basis, refCube);
		this->ReferenceElements.insert(std::make_pair(StandardElementCode::Cube, refCube));
	}

	bool IsGlobalBasis() { return false; }

	/*double VolumicTerm(Element* element, IBasisFunction3D* phi1, IBasisFunction3D* phi2)
	{
		Cube* cube = static_cast<Cube*>(element);
		return this->VolumicTerm(cube, phi1, phi2);
	}*/

	/*double VolumicTerm(Cube* element, IBasisFunction3D* phi1, IBasisFunction3D* phi2)
	{
		double h = element->Width;
		assert(h >= 0);
		DefInterval refInterval = phi1->DefinitionInterval();

		function<double(double, double, double)> functionToIntegrate = [phi1, phi2](double t, double u, double v) {
			return InnerProduct(phi1->Grad(t, u, v), phi2->Grad(t, u, v));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		double factor = h / refInterval.Length;
		return factor * Utils::Integral(nQuadPoints, functionToIntegrate, refInterval, refInterval, refInterval);
	}*/

	double CouplingTerm(Face* face, Element* element1, IBasisFunction3D* phi1, Element* element2, IBasisFunction3D* phi2)
	{
		assert(face->IsBetween(element1, element2));

		Cube* cube1 = static_cast<Cube*>(element1);
		Cube* cube2 = static_cast<Cube*>(element2);
		return this->CouplingTerm(face, cube1, phi1, cube2, phi2);
	}

	double CouplingTerm(Face* face, Cube* element1, IBasisFunction3D* phi1, Cube* element2, IBasisFunction3D* phi2)
	{
		DefInterval refInterval = phi1->DefinitionInterval();

		double h = element1->Width;

		auto n1 = element1->OuterNormalVector(face);
		auto n2 = element2->OuterNormalVector(face);

		double meanFactor = face->IsDomainBoundary ? 1 : 0.5;

		Face3D* face3D = (Face3D*)face;

		std::function<double(double, double)> functionToIntegrate;
		if (face3D->IsInXOYPlan)
		{
			double v1 = face == element1->TopFace ? refInterval.Right : refInterval.Left;
			double v2 = face == element2->TopFace ? refInterval.Right : refInterval.Left;

			functionToIntegrate = [n1, n2, phi1, phi2, v1, v2](double t, double u) {
				double meanGradPhi1_scal_jumpPhi2 = InnerProduct(phi1->Grad(t, u, v1), n2) * phi2->Eval(t, u, v2);
				double meanGradPhi2_scal_jumpPhi1 = InnerProduct(phi2->Grad(t, u, v2), n1) * phi1->Eval(t, u, v1);
				return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
			};
		}
		else if (face3D->IsInXOZPlan)
		{
			double u1 = face == element1->BackFace ? refInterval.Right : refInterval.Left;
			double u2 = face == element2->BackFace ? refInterval.Right : refInterval.Left;

			functionToIntegrate = [n1, n2, phi1, phi2, u1, u2](double t, double v) {
				double meanGradPhi1_scal_jumpPhi2 = InnerProduct(phi1->Grad(t, u1, v), n2) * phi2->Eval(t, u2, v);
				double meanGradPhi2_scal_jumpPhi1 = InnerProduct(phi2->Grad(t, u2, v), n1) * phi1->Eval(t, u1, v);
				return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
			};
		}
		else if (face3D->IsInYOZPlan)
		{
			double t1 = face == element1->RightFace ? refInterval.Right : refInterval.Left;
			double t2 = face == element2->RightFace ? refInterval.Right : refInterval.Left;

			functionToIntegrate = [n1, n2, phi1, phi2, t1, t2](double u, double v) {
				double meanGradPhi1_scal_jumpPhi2 = InnerProduct(phi1->Grad(t1, u, v), n2) * phi2->Eval(t2, u, v);
				double meanGradPhi2_scal_jumpPhi1 = InnerProduct(phi2->Grad(t2, u, v), n1) * phi1->Eval(t1, u, v);
				return meanGradPhi1_scal_jumpPhi2 + meanGradPhi2_scal_jumpPhi1;
			};
		}
		else
			return 0;

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 1;
		double factor = h / refInterval.Length;
		return -meanFactor * factor * Utils::Integral(nQuadPoints, functionToIntegrate, refInterval, refInterval);
	}

	double PenalizationTerm(Face* face, Element* element1, IBasisFunction3D* phi1, Element* element2, IBasisFunction3D* phi2, double penalizationCoefficient)
	{
		assert(face->IsBetween(element1, element2));

		Cube* cube1 = static_cast<Cube*>(element1);
		Cube* cube2 = static_cast<Cube*>(element2);
		return this->PenalizationTerm(face, cube1, phi1, cube2, phi2, penalizationCoefficient);
	}

	double PenalizationTerm(Face* face, Cube* element1, IBasisFunction3D* phi1, Cube* element2, IBasisFunction3D* phi2, double penalizationCoefficient)
	{
		double h = element1->Width;
		auto n1 = element1->OuterNormalVector(face);
		auto n2 = element2->OuterNormalVector(face);

		assert(InnerProduct(n1, n2) == 1 || InnerProduct(n1, n2) == -1);

		DefInterval refInterval = phi1->DefinitionInterval();

		Face3D* face3D = (Face3D*)face;

		function<double(double, double)> functionToIntegrate;

		if (face3D->IsInXOYPlan)
		{
			double v1 = face == element1->TopFace ? refInterval.Right : refInterval.Left;
			double v2 = face == element2->TopFace ? refInterval.Right : refInterval.Left;

			functionToIntegrate = [phi1, phi2, n1, n2, v1, v2](double t, double u) {
				return InnerProduct(n1, n2) * phi1->Eval(t, u, v1) * phi2->Eval(t, u, v2);
			};
		}
		else if (face3D->IsInXOZPlan)
		{
			double u1 = face == element1->BackFace ? refInterval.Right : refInterval.Left;
			double u2 = face == element2->BackFace ? refInterval.Right : refInterval.Left;

			functionToIntegrate = [phi1, phi2, n1, n2, u1, u2](double t, double v) {
				return InnerProduct(n1, n2) * phi1->Eval(t, u1, v) * phi2->Eval(t, u2, v);
			};
		}
		else if (face3D->IsInYOZPlan)
		{
			double t1 = face == element1->RightFace ? refInterval.Right : refInterval.Left;
			double t2 = face == element2->RightFace ? refInterval.Right : refInterval.Left;

			functionToIntegrate = [phi1, phi2, n1, n2, t1, t2](double u, double v) {
				return InnerProduct(n1, n2) * phi1->Eval(t1, u, v) * phi2->Eval(t2, u, v);
			};
		}
		else
			assert(false);

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		double jacobian = pow(h / refInterval.Length, 2);
		double integralJump1ScalarJump2 = jacobian * Utils::Integral(nQuadPoints, functionToIntegrate, refInterval, refInterval);
		return penalizationCoefficient * integralJump1ScalarJump2;
	}

	double RightHandSide(Element* element, IBasisFunction3D* phi)
	{
		Cube* cube = static_cast<Cube*>(element);
		return this->RightHandSide(cube, phi);
	}

	double RightHandSide(Cube* element, IBasisFunction3D* phi)
	{
		double x1 = element->X;
		double x2 = element->X + element->Width;
		double y1 = element->Y;
		double y2 = element->Y + element->Width;
		double z1 = element->Z;
		double z2 = element->Z + element->Width;

		DefInterval refInterval = phi->DefinitionInterval();

		function<double(double, double, double)> sourceTimesBasisFunction = NULL;
		if (refInterval.Left == -1 && refInterval.Right == 1)
		{
			sourceTimesBasisFunction = [this, phi, x1, x2, y1, y2, z1, z2](double t, double u, double v) {
				return this->_sourceFunction((x2 - x1) / 2 * t + (x2 + x1) / 2, (y2 - y1) / 2 * u + (y2 + y1) / 2, (z2 - z1) / 2 * v + (z2 + z1) / 2) * phi->Eval(t, u, v);
			};
		}
		else
		{
			sourceTimesBasisFunction = [this, phi, x1, x2, y1, y2, z1, z2](double t, double u, double v) {
				return this->_sourceFunction((x2 - x1) * t + x1, (y2 - y1) * u + y1, (z2 - z1) * v + z1) * phi->Eval(t, u, v);
			};
		}

		double jacobian = (x2 - x1) * (y2 - y1) * (z2 - z1) / pow(refInterval.Length, 3);
		return jacobian * Utils::Integral(sourceTimesBasisFunction, refInterval, refInterval, refInterval);
	}

private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2];
	}
};