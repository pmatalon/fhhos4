#pragma once
#include "CartesianElement.h"


class Cube : public CartesianElement<3>
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
		CartesianElement(number, DomPoint(x,y,z), width)
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
		cout << "Error: TODO use GradTransformation()." << endl;
		exit(EXIT_FAILURE);

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

};