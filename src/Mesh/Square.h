#pragma once
#include "Element.h"
#include "IntervalFace.h"
#include "../DG/Poisson_DG_Element.h"
#include "../DG/Poisson_DG_ReferenceElement.h"
#include "../HHO/Poisson_HHO_Element.h"
#include "../Utils/SourceFunction.h"
#include <assert.h>
using namespace std;

class Square : public Element<2>, public Poisson_DG_Element<2>, public Poisson_HHO_Element<2>
{
public:
	double X;
	double Y;
	double Width;

	Face<2>* NorthFace;
	Face<2>* SouthFace;
	Face<2>* EastFace;
	Face<2>* WestFace;

	Square(int number, double x, double y, double width) : Element(number)
	{
		this->X = x;
		this->Y = y;
		this->Width = width;
	}

	StandardElementCode StdElementCode()
	{
		return StandardElementCode::Square;
	}

	void SetNorthInterface(IntervalFace* face)
	{
		this->Faces.push_back(face);
		this->NorthFace = face;
	}

	void SetSouthInterface(IntervalFace* face)
	{
		this->Faces.push_back(face);
		this->SouthFace = face;
	}

	void SetEastInterface(IntervalFace* face)
	{
		this->Faces.push_back(face);
		this->EastFace = face;
	}

	void SetWestInterface(IntervalFace* face)
	{
		this->Faces.push_back(face);
		this->WestFace = face;
	}

	double* OuterNormalVector(Face<2>* face)
	{
		if (face == this->NorthFace)
			return new double[2]{ 0, 1 };
		if (face == this->SouthFace)
			return new double[2]{ 0, -1 };
		if (face == this->EastFace)
			return new double[2]{ 1, 0 };
		if (face == this->WestFace)
			return new double[2]{ -1, 0 };
		return NULL;
	}

	double DiffusionCoefficient(DiffusionPartition diffusionPartition)
	{
		return diffusionPartition.Coefficient(Point(this->X, this->Y));
	}

	double Integral(function<double(Point)> func)
	{
		double x1 = this->X;
		double x2 = this->X + this->Width;
		double y1 = this->Y;
		double y2 = this->Y + this->Width;

		return Utils::Integral(func, x1, x2, y1, y2);
	}

	double L2ErrorPow2(function<double(Point)> approximate, function<double(Point)> exactSolution)
	{
		double x1 = this->X;
		double x2 = this->X + this->Width;
		double y1 = this->Y;
		double y2 = this->Y + this->Width;

		function<double(double, double)> errorFunction = [exactSolution, approximate, x1, x2, y1, y2](double t, double u) {
			Point p;
			p.X = (x2 - x1) / 2 * t + (x2 + x1) / 2;
			p.Y = (y2 - y1) / 2 * u + (y2 + y1) / 2;
			return pow(exactSolution(p) - approximate(Point(t, u)), 2);
		};

		return (x2 - x1) * (y2 - y1) / 4 * Utils::Integral(errorFunction, -1, 1, -1, 1);
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double VolumicTerm(BasisFunction<2>* phi1, BasisFunction<2>* phi2, Poisson_DG_ReferenceElement<2>* referenceElement, DiffusionPartition diffusionPartition)
	{
		double kappa = this->DiffusionCoefficient(diffusionPartition);
		return kappa * referenceElement->VolumicTerm(phi1, phi2);
	}

	double MassTerm(BasisFunction<2>* phi1, BasisFunction<2>* phi2, Poisson_DG_ReferenceElement<2>* referenceElement)
	{
		double h = this->Width;
		return pow(h, 2) / 4 * referenceElement->MassTerm(phi1, phi2);
	}

	double SourceTerm(BasisFunction<2>* phi, SourceFunction* f)
	{
		double x1 = this->X;
		double x2 = this->X + this->Width;
		double y1 = this->Y;
		double y2 = this->Y + this->Width;

		function<double(double, double)> sourceTimesBasisFunction = [f, phi, x1, x2, y1, y2](double t, double u) {
			Point p;
			p.X = (x2 - x1) / 2 * t + (x2 + x1) / 2;
			p.Y = (y2 - y1) / 2 * u + (y2 + y1) / 2;
			return f->Eval(p) * phi->Eval(Point(t, u));
		};

		double jacobian = (x2 - x1) * (y2 - y1) / 4;
		return jacobian * Utils::Integral(sourceTimesBasisFunction, -1,1, -1,1);
	}

	function<double(Point)> EvalPhiOnFace(Face<2>* face, BasisFunction<2>* p_phi)
	{
		IBasisFunction2D* phi = static_cast<IBasisFunction2D*>(p_phi);

		function<double(Point)> evalOnFace = NULL;
		if (face == this->EastFace || face == this->WestFace)
		{
			double tFixed = face == this->EastFace ? 1 : -1;
			evalOnFace = [phi, tFixed](Point point1D) {
				double u = point1D.X;
				return phi->Eval(tFixed, u);
			};
		}
		else if (face == this->SouthFace || face == this->NorthFace)
		{
			double uFixed = face == this->NorthFace ? 1 : -1;
			evalOnFace = [phi, uFixed](Point point1D) {
				double t = point1D.X;
				return phi->Eval(t, uFixed);
			};
		}
		else
			assert(false);
		return evalOnFace;
	}


	function<double*(Point)> GradPhiOnFace(Face<2>* face, BasisFunction<2>* p_phi)
	{
		IBasisFunction2D* phi = static_cast<IBasisFunction2D*>(p_phi);

		function<double*(Point)> gradOnFace = NULL;
		if (face == this->EastFace || face == this->WestFace)
		{
			double tFixed = face == this->EastFace ? 1 : -1;
			gradOnFace = [phi, tFixed](Point point1D) {
				double u = point1D.X;
				return phi->Grad(tFixed, u);
			};
		}
		else if (face == this->SouthFace || face == this->NorthFace)
		{
			double uFixed = face == this->NorthFace ? 1 : -1;
			gradOnFace = [phi, uFixed](Point point1D) {
				double t = point1D.X;
				return phi->Grad(t, uFixed);
			};
		}
		else
			assert(false);
		return gradOnFace;
	}

	//-------------------------------------------------------------------//
	//                 Poisson_HHO_Element implementation                //
	//-------------------------------------------------------------------//

	double IntegralGradGrad(BasisFunction<2>* p_phi1, BasisFunction<2>* p_phi2)
	{
		IBasisFunction2D* phi1 = static_cast<IBasisFunction2D*>(p_phi1);
		IBasisFunction2D* phi2 = static_cast<IBasisFunction2D*>(p_phi2);

		function<double(double, double)> functionToIntegrate = [phi1, phi2](double t, double u) {
			return InnerProduct(phi1->Grad(t, u), phi2->Grad(t, u));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		return Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1, -1, 1);
	}

	double S(BasisFunction<2>* phiReconstruct1, BasisFunction<2>* phiReconstruct2)
	{
		return this->IntegralGradGrad(phiReconstruct1, phiReconstruct2);
	}

	double Bt(BasisFunction<2>* phiReconstruct, BasisFunction<2>* phiElement)
	{
		double integralGradGrad = this->IntegralGradGrad(phiReconstruct, phiElement);
		
		double sumFaces = 0;
		for (auto face : this->Faces)
		{
			auto phi = this->EvalPhiOnFace(face, phiElement);
			auto gradPhi = this->GradPhiOnFace(face, phiReconstruct);
			auto normal = this->OuterNormalVector(face);

			std::function<double(double)> functionToIntegrate = [phi, gradPhi, normal](double u) {
				Point p(u);
				return InnerProduct(gradPhi(p), normal) * phi(p);
			};

			int nQuadPoints = phiReconstruct->GetDegree() + phiElement->GetDegree() + 1;
			sumFaces += Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1);
		}
		
		return integralGradGrad - sumFaces;
	}

	double Bf(BasisFunction<2>* phiReconstruct, BasisFunction<1>* phiFace)
	{
		double sumFaces = 0;
		for (auto face : this->Faces)
		{
			auto gradPhi = this->GradPhiOnFace(face, phiReconstruct);
			auto normal = this->OuterNormalVector(face);

			std::function<double(double)> functionToIntegrate = [phiFace, gradPhi, normal](double u) {
				Point p(u);
				return InnerProduct(gradPhi(p), normal) * phiFace->Eval(p);
			};

			int nQuadPoints = phiReconstruct->GetDegree() + phiFace->GetDegree() + 1;
			sumFaces += Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1);
		}
		return sumFaces;
	}

/*private:
	static double InnerProduct(double* vector1, double* vector2)
	{
		return vector1[0] * vector2[0] + vector1[1] * vector2[1];
	}*/
};

//Square* Square::ReferenceSquare = new Square(-1, -1, -1, 2);