#pragma once
#include "Element.h"


class Cube : public Element, public Poisson_DG_Element<3>
{
public:
	double X;
	double Y;
	double Z;
	double Width;

	Face* TopFace;
	Face* BottomFace;
	Face* FrontFace;
	Face* BackFace;
	Face* LeftFace;
	Face* RightFace;
private:
	Cube* _neighbours;
public:
	Cube(int number, double x, double y, double z, double width) : Element(number)
	{
		this->X = x;
		this->Y = y;
		this->Z = z;
		this->Width = width;
	}

	StandardElementCode StdElementCode()
	{
		return StandardElementCode::Cube;
	}

	void SetTopInterface(Face* interface)
	{
		this->Faces.push_back(interface);
		this->TopFace = interface;
	}

	void SetBottomInterface(Face* interface)
	{
		this->Faces.push_back(interface);
		this->BottomFace = interface;
	}

	void SetFrontInterface(Face* interface)
	{
		this->Faces.push_back(interface);
		this->FrontFace = interface;
	}

	void SetBackInterface(Face* interface)
	{
		this->Faces.push_back(interface);
		this->BackFace = interface;
	}

	void SetLeftInterface(Face* interface)
	{
		this->Faces.push_back(interface);
		this->LeftFace = interface;
	}

	void SetRightInterface(Face* interface)
	{
		this->Faces.push_back(interface);
		this->RightFace = interface;
	}

	double* OuterNormalVector(Face* interface)
	{
		if (interface == this->TopFace)
			return new double[3]{ 0, 0, 1 };
		if (interface == this->BottomFace)
			return new double[3]{ 0, 0, -1 };
		if (interface == this->FrontFace)
			return new double[3]{ 0, -1, 0 };
		if (interface == this->BackFace)
			return new double[3]{ 0, 1, 0 };
		if (interface == this->LeftFace)
			return new double[3]{ -1, 0, 0 };
		if (interface == this->RightFace)
			return new double[3]{ 1, 0, 0 };
		return NULL;
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

	double L2ErrorPow2(function<double(Point)> approximate, function<double(Point)> exactSolution)
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
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double VolumicTerm(BasisFunction<3>* phi1, BasisFunction<3>* phi2, Poisson_DG_ReferenceElement<3>* referenceElement)
	{
		double h = this->Width;
		return h / 2 * referenceElement->VolumicTerm(phi1, phi2);
	}

	double MassTerm(BasisFunction<3>* phi1, BasisFunction<3>* phi2, Poisson_DG_ReferenceElement<3>* referenceElement)
	{
		double h = this->Width;
		return pow(h, 3) / 8 * referenceElement->MassTerm(phi1, phi2);
	}

	double SourceTerm(BasisFunction<3>* phi, SourceFunction* f)
	{
		double x1 = this->X;
		double x2 = this->X + this->Width;
		double y1 = this->Y;
		double y2 = this->Y + this->Width;
		double z1 = this->Z;
		double z2 = this->Z + this->Width;

		function<double(double, double, double)> sourceTimesBasisFunction = [f, phi, x1, x2, y1, y2, z1, z2](double t, double u, double v) {
			Point p;
			p.X = (x2 - x1) / 2 * t + (x2 + x1) / 2;
			p.Y = (y2 - y1) / 2 * u + (y2 + y1) / 2;
			p.Z = (z2 - z1) / 2 * v + (z2 + z1) / 2;
			return f->Eval(p) * phi->Eval(Point(t, u, v));
		};

		double jacobian = (x2 - x1) * (y2 - y1) * (z2 - z1) / 8;
		return jacobian * Utils::Integral(sourceTimesBasisFunction, -1,1, -1,1, -1,1);
	}

	function<double(Point)> EvalPhiOnFace(Face* face, BasisFunction<3>* p_phi)
	{
		IBasisFunction3D* phi = dynamic_cast<IBasisFunction3D*>(p_phi);

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


	function<double*(Point)> GradPhiOnFace(Face* face, BasisFunction<3>* p_phi)
	{
		IBasisFunction3D* phi = dynamic_cast<IBasisFunction3D*>(p_phi);

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
};