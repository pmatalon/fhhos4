#pragma once
#include "Element.h"


class Cube : public Element, public Poisson_DG_Element
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

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double VolumicTerm(BasisFunction* phi1, BasisFunction* phi2, Poisson_DG_ReferenceElement* referenceElement)
	{
		double h = this->Width;
		double factor = h / 2;
		return factor * referenceElement->VolumicTerm(phi1, phi2);
	}

	function<double(Point)> EvalPhiOnFace(Face* face, BasisFunction* p_phi)
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


	function<double*(Point)> GradPhiOnFace(Face* face, BasisFunction* p_phi)
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