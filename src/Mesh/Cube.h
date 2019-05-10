#pragma once
#include "CartesianElement.h"


class Cube : public CartesianElement<3>, public Poisson_DG_Element<3>, public Poisson_HHO_Element<3>
{
public:
	Face<3>* TopFace;
	Face<3>* BottomFace;
	Face<3>* FrontFace;
	Face<3>* BackFace;
	Face<3>* LeftFace;
	Face<3>* RightFace;

public:
	Cube(int number, double x, double y, double z, double width) : 
		Element(number), 
		CartesianElement(number, DomPoint(x,y,z), width), 
		Poisson_DG_Element(number), 
		Poisson_HHO_Element(number)
	{ }

	void SetTopInterface(Face<3>* face)
	{
		this->AddFace(face);
		this->TopFace = face;
	}

	void SetBottomInterface(Face<3>* face)
	{
		this->AddFace(face);
		this->BottomFace = face;
	}

	void SetFrontInterface(Face<3>* face)
	{
		this->AddFace(face);
		this->FrontFace = face;
	}

	void SetBackInterface(Face<3>* face)
	{
		this->AddFace(face);
		this->BackFace = face;
	}

	void SetLeftInterface(Face<3>* face)
	{
		this->AddFace(face);
		this->LeftFace = face;
	}

	void SetRightInterface(Face<3>* face)
	{
		this->AddFace(face);
		this->RightFace = face;
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	StandardElementCode StdElementCode()
	{
		return StandardElementCode::Cube;
	}

	vector<double> OuterNormalVector(Face<3>* face)
	{
		if (face == this->TopFace)
			return vector<double>{ 0, 0, 1 };
		if (face == this->BottomFace)
			return vector<double>{ 0, 0, -1 };
		if (face == this->FrontFace)
			return vector<double>{ 0, -1, 0 };
		if (face == this->BackFace)
			return vector<double>{ 0, 1, 0 };
		if (face == this->LeftFace)
			return vector<double>{ -1, 0, 0 };
		if (face == this->RightFace)
			return vector<double>{ 1, 0, 0 };
		assert(false);
	}

	double IntegralGlobalFunction(function<double(DomPoint)> func)
	{
		double x1 = this->Origin.X;
		double x2 = this->Origin.X + this->Width;
		double y1 = this->Origin.Y;
		double y2 = this->Origin.Y + this->Width;
		double z1 = this->Origin.Z;
		double z2 = this->Origin.Z + this->Width;

		return Utils::Integral(func, x1, x2, y1, y2, z1, z2);
	}

	function<double(Point)> EvalPhiOnFace(Face<3>* face, BasisFunction<3>* p_phi) override
	{
		IBasisFunction3D* phi = static_cast<IBasisFunction3D*>(p_phi);

		function<double(RefPoint)> evalOnFace = NULL;
		if (face == this->TopFace || face == this->BottomFace) // XOY plan
		{
			double vFixed = face == this->TopFace ? 1 : -1;
			evalOnFace = [phi, vFixed](RefPoint point2D) {
				double t = point2D.X;
				double u = point2D.Y;
				return phi->Eval(t, u, vFixed);
			};
		}
		else if (face == this->BackFace || face == this->FrontFace) // XOZ plan
		{
			double uFixed = face == this->BackFace ? 1 : -1;
			evalOnFace = [phi, uFixed](RefPoint point2D) {
				double t = point2D.X;
				double v = point2D.Y;
				return phi->Eval(t, uFixed, v);
			};
		}
		else if (face == this->RightFace || face == this->LeftFace) // YOZ plan
		{
			double tFixed = face == this->RightFace ? 1 : -1;
			evalOnFace = [phi, tFixed](RefPoint point2D) {
				double u = point2D.X;
				double v = point2D.Y;
				return phi->Eval(tFixed, u, v);
			};
		}
		else
			assert(false);
		return evalOnFace;
	}

	function<vector<double>(RefPoint)> GradPhiOnFace(Face<3>* face, BasisFunction<3>* p_phi) override
	{
		IBasisFunction3D* phi = static_cast<IBasisFunction3D*>(p_phi);

		function<vector<double>(RefPoint)> gradOnFace = NULL;
		if (face == this->TopFace || face == this->BottomFace) // XOY plan
		{
			double vFixed = face == this->TopFace ? 1 : -1;
			gradOnFace = [phi, vFixed](RefPoint point2D) {
				double t = point2D.X;
				double u = point2D.Y;
				return phi->Grad(t, u, vFixed);
			};
		}
		else if (face == this->BackFace || face == this->FrontFace) // XOZ plan
		{
			double uFixed = face == this->BackFace ? 1 : -1;
			gradOnFace = [phi, uFixed](RefPoint point2D) {
				double t = point2D.X;
				double v = point2D.Y;
				return phi->Grad(t, uFixed, v);
			};
		}
		else if (face == this->RightFace || face == this->LeftFace) // YOZ plan
		{
			double tFixed = face == this->RightFace ? 1 : -1;
			gradOnFace = [phi, tFixed](RefPoint point2D) {
				double u = point2D.X;
				double v = point2D.Y;
				return phi->Grad(tFixed, u, v);
			};
		}
		else
			assert(false);
		return gradOnFace;
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double VolumicTerm(BasisFunction<3>* phi1, BasisFunction<3>* phi2, Poisson_DG_ReferenceElement<3>* referenceElement, DiffusionPartition diffusionPartition)
	{
		double h = this->Width;
		double kappa = CartesianElement::DiffusionCoefficient(diffusionPartition);
		return h / 2 * kappa * referenceElement->VolumicTerm(phi1, phi2);
	}

	double MassTerm(BasisFunction<3>* phi1, BasisFunction<3>* phi2, Poisson_DG_ReferenceElement<3>* referenceElement)
	{
		double h = this->Width;
		return pow(h, 3) / 8 * referenceElement->MassTerm(phi1, phi2);
	}

	double SourceTerm(BasisFunction<3>* phi, SourceFunction* f)
	{
		return CartesianElement::SourceTerm(phi, f);
	}

	//-------------------------------------------------------------------//
	//                 Poisson_HHO_Element implementation                //
	//-------------------------------------------------------------------//

	double St(BasisFunction<3>* reconstructPhi1, BasisFunction<3>* reconstructPhi2)
	{
		return this->IntegralGradGrad(reconstructPhi1, reconstructPhi2);
	}

	double Lt(BasisFunction<3>* phi)
	{
		return CartesianElement::Integral(phi);
	}

	double Bt(BasisFunction<3>* reconstructPhi, BasisFunction<3>* cellPhi)
	{
		if (reconstructPhi->GetDegree() == 0)
			return 0;


		double integralGradGrad = this->IntegralGradGrad(reconstructPhi, cellPhi);
		double h = this->Width;

		double sumFaces = 0;
		for (auto face : this->Faces)
		{
			auto phi = this->EvalPhiOnFace(face, cellPhi);
			auto gradPhi = this->GradPhiOnFace(face, reconstructPhi);
			auto normal = this->OuterNormalVector(face);

			std::function<double(double, double)> functionToIntegrate = [phi, gradPhi, normal](double u, double v) {
				RefPoint p(u, v);
				return InnerProduct(gradPhi(p), normal) * phi(p);
			};

			int nQuadPoints = reconstructPhi->GetDegree() + cellPhi->GetDegree() + 1;
			double integralFace = h/2 * Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1, -1, 1);
			sumFaces += integralFace;
		}

		return integralGradGrad - sumFaces;
	}

	double Bf(BasisFunction<3>* reconstructPhi, BasisFunction<2>* facePhi, Face<3>* face)
	{
		if (reconstructPhi->GetDegree() == 0)
			return 0;

		auto gradPhi = this->GradPhiOnFace(face, reconstructPhi);
		auto normal = this->OuterNormalVector(face);

		std::function<double(double, double)> functionToIntegrate = [facePhi, gradPhi, normal](double u, double v) {
			RefPoint p(u, v);
			return InnerProduct(gradPhi(p), normal) * facePhi->Eval(p);
		};

		double h = this->Width;
		int nQuadPoints = reconstructPhi->GetDegree() + facePhi->GetDegree() + 1;
		return h / 2 * Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1, -1, 1);
	}
};