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

	DomPoint BottomLeftCorner;
	DomPoint TopLeftCorner;
	DomPoint TopRightCorner;
	DomPoint BottomRightCorner;

	Square(int number, double x, double y, double width) : 
		Element(number), 
		CartesianElement(number, DomPoint(x, y), width), 
		Poisson_DG_Element(number), Poisson_HHO_Element(number),
		BottomLeftCorner(x, y),
		TopLeftCorner(x, y + width),
		TopRightCorner(x + width, y + width),
		BottomRightCorner(x + width, y)
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

	double IntegralGlobalFunction(function<double(DomPoint)> func) override
	{
		double x1 = this->Origin.X;
		double x2 = this->Origin.X + this->Width;
		double y1 = this->Origin.Y;
		double y2 = this->Origin.Y + this->Width;

		return Utils::Integral(func, x1, x2, y1, y2);
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
				RefPoint p(u);
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
			RefPoint p(u);
			return InnerProduct(gradPhi(p), normal) * facePhi->Eval(p);
		};

		int nQuadPoints = reconstructPhi->GetDegree() + facePhi->GetDegree() + 1;
		return Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1);
	}
};