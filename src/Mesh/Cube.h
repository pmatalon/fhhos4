#pragma once
#include "CartesianElement.h"


class Cube : public CartesianElement<3>, public Poisson_DG_Element<3>, public Poisson_HHO_Element<3>
{
public:
	double X;
	double Y;
	double Z;
	double Width;

	Face<3>* TopFace;
	Face<3>* BottomFace;
	Face<3>* FrontFace;
	Face<3>* BackFace;
	Face<3>* LeftFace;
	Face<3>* RightFace;

	Reconstructor<3>* HHOReconstructor = NULL;

public:
	Cube(int number, double x, double y, double z, double width) : CartesianElement(number, width)
	{
		this->X = x;
		this->Y = y;
		this->Z = z;
		this->Width = width;
	}

	/*double GetDiameter()
	{
		return this->Width;
	}*/

	StandardElementCode StdElementCode()
	{
		return StandardElementCode::Cube;
	}

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

	double* OuterNormalVector(Face<3>* face)
	{
		if (face == this->TopFace)
			return new double[3]{ 0, 0, 1 };
		if (face == this->BottomFace)
			return new double[3]{ 0, 0, -1 };
		if (face == this->FrontFace)
			return new double[3]{ 0, -1, 0 };
		if (face == this->BackFace)
			return new double[3]{ 0, 1, 0 };
		if (face == this->LeftFace)
			return new double[3]{ -1, 0, 0 };
		if (face == this->RightFace)
			return new double[3]{ 1, 0, 0 };
		return NULL;
	}

	double DiffusionCoefficient(DiffusionPartition diffusionPartition)
	{
		return diffusionPartition.Coefficient(Point(this->X, this->Y, this->Z));
	}

	double Integral(function<double(Point)> func)
	{
		double x1 = this->X;
		double x2 = this->X + this->Width;
		double y1 = this->Y;
		double y2 = this->Y + this->Width;
		double z1 = this->Z;
		double z2 = this->Z + this->Width;

		return Utils::Integral(func, x1, x2, y1, y2, z1, z2);
	}

	/*double Integral(BasisFunction<3>* phi)
	{
		double x1 = this->X;
		double x2 = this->X + this->Width;
		double y1 = this->Y;
		double y2 = this->Y + this->Width;
		double z1 = this->Z;
		double z2 = this->Z + this->Width;

		return (x2 - x1) * (y2 - y1) * (z2 - z1) / 8 * Utils::Integral(phi, -1, 1, -1, 1, -1, 1);
	}*/

	Point ConvertToDomain(Point referenceElementPoint) override
	{
		double x1 = this->X;
		double x2 = this->X + this->Width;
		double y1 = this->Y;
		double y2 = this->Y + this->Width;
		double z1 = this->Z;
		double z2 = this->Z + this->Width;

		double t = referenceElementPoint.X;
		double u = referenceElementPoint.Y;
		double v = referenceElementPoint.Z;

		Point p;
		p.X = (x2 - x1) / 2 * t + (x2 + x1) / 2;
		p.Y = (y2 - y1) / 2 * u + (y2 + y1) / 2;
		p.Z = (z2 - z1) / 2 * v + (z2 + z1) / 2;
		return p;
	}

	/*double L2ErrorPow2(function<double(Point)> approximate, function<double(Point)> exactSolution) override
	{
		double x1 = this->X;
		double x2 = this->X + this->Width;
		double y1 = this->Y;
		double y2 = this->Y + this->Width;
		double z1 = this->Z;
		double z2 = this->Z + this->Width;

		function<double(double, double, double)> errorFunction = [exactSolution, approximate, x1, x2, y1, y2, z1, z2](double t, double u, double v) {
			Point p;
			p.X = (x2 - x1) / 2 * t + (x2 + x1) / 2;
			p.Y = (y2 - y1) / 2 * u + (y2 + y1) / 2;
			p.Z = (z2 - z1) / 2 * v + (z2 + z1) / 2;
			return pow(exactSolution(p) - approximate(Point(t, u, v)), 2);
		};

		return (x2 - x1) * (y2 - y1) * (z2 - z1) / 8 * Utils::Integral(errorFunction, -1, 1, -1, 1, -1, 1);
	}*/

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double VolumicTerm(BasisFunction<3>* phi1, BasisFunction<3>* phi2, Poisson_DG_ReferenceElement<3>* referenceElement, DiffusionPartition diffusionPartition)
	{
		double h = this->Width;
		double kappa = this->DiffusionCoefficient(diffusionPartition);
		return h / 2 * kappa * referenceElement->VolumicTerm(phi1, phi2);
	}

	double MassTerm(BasisFunction<3>* phi1, BasisFunction<3>* phi2, Poisson_DG_ReferenceElement<3>* referenceElement)
	{
		double h = this->Width;
		return pow(h, 3) / 8 * referenceElement->MassTerm(phi1, phi2);
	}

	double MassTerm(BasisFunction<3>* p_phi1, BasisFunction<3>* p_phi2) override
	{
		/*double h = this->Width;

		IBasisFunction3D* phi1 = static_cast<IBasisFunction3D*>(p_phi1);
		IBasisFunction3D* phi2 = static_cast<IBasisFunction3D*>(p_phi2);

		function<double(double, double, double)> functionToIntegrate = [phi1, phi2](double t, double u, double v) {
			return phi1->Eval(t, u, v)*phi2->Eval(t, u, v);
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		return pow(h, 3) / 8 * Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1, -1, 1, -1, 1);*/
		return CartesianElement::MassTerm(p_phi1, p_phi2);
	}
	
	Eigen::MatrixXd MassMatrix(FunctionalBasis<3>* basis)
	{
		return Element::MassMatrix(basis);
	}

	Eigen::MatrixXd MassMatrix(FunctionalBasis<3>* basis1, FunctionalBasis<3>* basis2)
	{
		return Element::MassMatrix(basis1, basis2);
	}

	double SourceTerm(BasisFunction<3>* phi, SourceFunction* f)
	{
		return CartesianElement::SourceTerm(phi, f);
	}

	function<double(Point)> EvalPhiOnFace(Face<3>* face, BasisFunction<3>* p_phi)
	{
		IBasisFunction3D* phi = static_cast<IBasisFunction3D*>(p_phi);

		function<double(Point)> evalOnFace = NULL;
		if (face == this->TopFace || face == this->BottomFace) // XOY plan
		{
			double vFixed = face == this->TopFace ? 1 : -1;
			evalOnFace = [phi, vFixed](Point point2D) {
				double t = point2D.X;
				double u = point2D.Y;
				return phi->Eval(t, u, vFixed);
			};
		}
		else if (face == this->BackFace || face == this->FrontFace) // XOZ plan
		{
			double uFixed = face == this->BackFace ? 1 : -1;
			evalOnFace = [phi, uFixed](Point point2D) {
				double t = point2D.X;
				double v = point2D.Y;
				return phi->Eval(t, uFixed, v);
			};
		}
		else if (face == this->RightFace || face == this->LeftFace) // YOZ plan
		{
			double tFixed = face == this->RightFace ? 1 : -1;
			evalOnFace = [phi, tFixed](Point point2D) {
				double u = point2D.X;
				double v = point2D.Y;
				return phi->Eval(tFixed, u, v);
			};
		}
		else
			assert(false);
		return evalOnFace;
	}


	function<double*(Point)> GradPhiOnFace(Face<3>* face, BasisFunction<3>* p_phi)
	{
		IBasisFunction3D* phi = static_cast<IBasisFunction3D*>(p_phi);

		function<double*(Point)> gradOnFace = NULL;
		if (face == this->TopFace || face == this->BottomFace) // XOY plan
		{
			double vFixed = face == this->TopFace ? 1 : -1;
			gradOnFace = [phi, vFixed](Point point2D) {
				double t = point2D.X;
				double u = point2D.Y;
				return phi->Grad(t, u, vFixed);
			};
		}
		else if (face == this->BackFace || face == this->FrontFace) // XOZ plan
		{
			double uFixed = face == this->BackFace ? 1 : -1;
			gradOnFace = [phi, uFixed](Point point2D) {
				double t = point2D.X;
				double v = point2D.Y;
				return phi->Grad(t, uFixed, v);
			};
		}
		else if (face == this->RightFace || face == this->LeftFace) // YOZ plan
		{
			double tFixed = face == this->RightFace ? 1 : -1;
			gradOnFace = [phi, tFixed](Point point2D) {
				double u = point2D.X;
				double v = point2D.Y;
				return phi->Grad(tFixed, u, v);
			};
		}
		else
			assert(false);
		return gradOnFace;
	}

	//-------------------------------------------------------------------//
	//                 Poisson_HHO_Element implementation                //
	//-------------------------------------------------------------------//

	void InitReconstructor(FunctionalBasis<3>* reconstructionBasis, FunctionalBasis<3>* elementBasis, FunctionalBasis<2>* faceBasis)
	{
		this->HHOReconstructor = new Reconstructor<3>(this, reconstructionBasis, elementBasis, faceBasis);
	}

	Eigen::VectorXd Reconstruct(Eigen::VectorXd hybridVector)
	{
		return this->HHOReconstructor->Reconstruct(hybridVector);
	}

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
				Point p(u, v);
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
			Point p(u, v);
			return InnerProduct(gradPhi(p), normal) * facePhi->Eval(p);
		};

		double h = this->Width;
		int nQuadPoints = reconstructPhi->GetDegree() + facePhi->GetDegree() + 1;
		return h / 2 * Utils::Integral(nQuadPoints, functionToIntegrate, -1, 1, -1, 1);
	}

	int LocalNumberOf(Face<3>* face)
	{
		return Element::LocalNumberOf(face);
	}

	int FirstDOFLocalNumber(Face<3>* face)
	{
		return this->HHOReconstructor->FirstDOFNumber(face);
	}

	double ConsistencyTerm(BasisFunction<3>* cellPhi1, BasisFunction<3>* cellPhi2)
	{
		return this->HHOReconstructor->ConsistencyTerm(cellPhi1, cellPhi2);
	}
	double ConsistencyTerm(Face<3>* face, BasisFunction<3>* cellPhi, BasisFunction<2>* facePhi)
	{
		return this->HHOReconstructor->ConsistencyTerm(face, cellPhi, facePhi);
	}
	double ConsistencyTerm(Face<3>* face1, BasisFunction<2>* facePhi1, Face<3>* face2, BasisFunction<2>* facePhi2)
	{
		return this->HHOReconstructor->ConsistencyTerm(face1, facePhi1, face2, facePhi2);
	}

	double StabilizationTerm(BasisFunction<3>* cellPhi1, BasisFunction<3>* cellPhi2)
	{
		return this->HHOReconstructor->StabilizationTerm(cellPhi1, cellPhi2);
	}
	double StabilizationTerm(Face<3>* face, BasisFunction<3>* cellPhi, BasisFunction<2>* facePhi)
	{
		return this->HHOReconstructor->StabilizationTerm(face, cellPhi, facePhi);
	}
	double StabilizationTerm(Face<3>* face1, BasisFunction<2>* facePhi1, Face<3>* face2, BasisFunction<2>* facePhi2)
	{
		return this->HHOReconstructor->StabilizationTerm(face1, facePhi1, face2, facePhi2);
	}

	virtual double ReconstructionTerm(BasisFunction<3>* reconstrucPhi, BasisFunction<3>* cellPhi)
	{
		return this->HHOReconstructor->ReconstructionTerm(reconstrucPhi, cellPhi);
	}
	virtual double ReconstructionTerm(BasisFunction<3>* reconstrucPhi, Face<3>* face, BasisFunction<2>* facePhi)
	{
		return this->HHOReconstructor->ReconstructionTerm(reconstrucPhi, face, facePhi);
	}
};