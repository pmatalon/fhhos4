#pragma once
#include "CartesianElement.h"
#include "IntervalFace.h"
#include "../DG/Poisson_DG_Element.h"
#include "../DG/Poisson_DG_ReferenceElement.h"
#include "../HHO/Poisson_HHO_Element.h"
#include "../Utils/SourceFunction.h"
#include <assert.h>
using namespace std;

class Square : public CartesianElement<2>, public Poisson_DG_Element<2>, public Poisson_HHO_Element<2>
{
public:
	Face<2>* NorthFace;
	Face<2>* SouthFace;
	Face<2>* EastFace;
	Face<2>* WestFace;

	Square(int number, double x, double y, double width) : Element(number), CartesianElement(number, Point(x, y), width), Poisson_DG_Element(number), Poisson_HHO_Element(number)
	{
	}

	void SetNorthInterface(IntervalFace* face)
	{
		this->AddFace(face);
		this->NorthFace = face;
	}

	void SetSouthInterface(IntervalFace* face)
	{
		this->AddFace(face);
		this->SouthFace = face;
	}

	void SetEastInterface(IntervalFace* face)
	{
		this->AddFace(face);
		this->EastFace = face;
	}

	void SetWestInterface(IntervalFace* face)
	{
		this->AddFace(face);
		this->WestFace = face;
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	StandardElementCode StdElementCode()
	{
		return StandardElementCode::Square;
	}

	vector<double> OuterNormalVector(Face<2>* face)
	{
		if (face == this->NorthFace)
			return vector<double>{ 0, 1 };
		if (face == this->SouthFace)
			return vector<double>{ 0, -1 };
		if (face == this->EastFace)
			return vector<double>{ 1, 0 };
		if (face == this->WestFace)
			return vector<double>{ -1, 0 };
		assert(false);
	}

	double IntegralGlobalFunction(function<double(Point)> func) override
	{
		double x1 = this->Origin.X;
		double x2 = this->Origin.X + this->Width;
		double y1 = this->Origin.Y;
		double y2 = this->Origin.Y + this->Width;

		return Utils::Integral(func, x1, x2, y1, y2);
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

	function<vector<double>(Point)> GradPhiOnFace(Face<2>* face, BasisFunction<2>* p_phi)
	{
		IBasisFunction2D* phi = static_cast<IBasisFunction2D*>(p_phi);

		function<vector<double>(Point)> gradOnFace = NULL;
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

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double VolumicTerm(BasisFunction<2>* phi1, BasisFunction<2>* phi2, Poisson_DG_ReferenceElement<2>* referenceElement, DiffusionPartition diffusionPartition)
	{
		double kappa = CartesianElement::DiffusionCoefficient(diffusionPartition);
		return kappa * referenceElement->VolumicTerm(phi1, phi2);
	}

	double MassTerm(BasisFunction<2>* phi1, BasisFunction<2>* phi2, Poisson_DG_ReferenceElement<2>* referenceElement)
	{
		double h = this->Width;
		return pow(h, 2) / 4 * referenceElement->MassTerm(phi1, phi2);
	}

	double SourceTerm(BasisFunction<2>* phi, SourceFunction* f)
	{
		return CartesianElement::SourceTerm(phi, f);
	}

	//-------------------------------------------------------------------//
	//                 Poisson_HHO_Element implementation                //
	//-------------------------------------------------------------------//
	
	double St(BasisFunction<2>* reconstructPhi1, BasisFunction<2>* reconstructPhi2)
	{
		return this->IntegralGradGrad(reconstructPhi1, reconstructPhi2);
	}

	double Lt(BasisFunction<2>* phi)
	{
		return CartesianShape::Integral(phi);
	}

	double Bt(BasisFunction<2>* reconstructPhi, BasisFunction<2>* cellPhi)
	{
		if (reconstructPhi->GetDegree() == 0)
			return 0;
		
		double integralGradGrad = this->IntegralGradGrad(reconstructPhi, cellPhi);
		
		double sumFaces = 0;
		for (auto face : this->Faces)
		{
			auto phi = this->EvalPhiOnFace(face, cellPhi);
			auto gradPhi = this->GradPhiOnFace(face, reconstructPhi);
			auto normal = this->OuterNormalVector(face);

			std::function<double(double)> functionToIntegrate = [phi, gradPhi, normal](double u) {
				Point p(u);
				return InnerProduct(gradPhi(p), normal) * phi(p);
			};

			int nQuadPoints = reconstructPhi->GetDegree() + cellPhi->GetDegree() + 1;
			double integralFace = Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1);
			sumFaces += integralFace;
		}
		
		return integralGradGrad - sumFaces;
	}

	double Bf(BasisFunction<2>* reconstructPhi, BasisFunction<1>* facePhi, Face<2>* face)
	{
		if (reconstructPhi->GetDegree() == 0)
			return 0;

		auto gradPhi = this->GradPhiOnFace(face, reconstructPhi);
		auto normal = this->OuterNormalVector(face);

		std::function<double(double)> functionToIntegrate = [facePhi, gradPhi, normal](double u) {
			Point p(u);
			return InnerProduct(gradPhi(p), normal) * facePhi->Eval(p);
		};

		int nQuadPoints = reconstructPhi->GetDegree() + facePhi->GetDegree() + 1;
		return Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1);
	}
};