#pragma once
#include "Element.h"
#include "Face2D.h"
#include "Poisson_DG_Element.h"
#include "Poisson_DG_ReferenceElement.h"
#include "SourceFunction.h"
#include <assert.h>

class Square : public Element, public Poisson_DG_Element<2>
{
public:
	double X;
	double Y;
	double Width;

	Face2D* NorthFace;
	Face2D* SouthFace;
	Face2D* EastFace;
	Face2D* WestFace;
private:
	//Square* _neighbours;
public:
	//static Square* ReferenceSquare;

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

	void SetNorthInterface(Face2D* face)
	{
		this->Faces.push_back(face);
		this->NorthFace = face;
		face->X1 = this->X;
		face->Y1 = this->Y + this->Width;
		face->X2 = this->X + this->Width;
		face->Y2 = this->Y + this->Width;
	}

	void SetSouthInterface(Face2D* face)
	{
		this->Faces.push_back(face);
		this->SouthFace = face;
		face->X1 = this->X;
		face->Y1 = this->Y;
		face->X2 = this->X + this->Width;
		face->Y2 = this->Y;
	}

	void SetEastInterface(Face2D* face)
	{
		this->Faces.push_back(face);
		this->EastFace = face;
		face->X1 = this->X + this->Width;
		face->Y1 = this->Y;
		face->X2 = this->X + this->Width;
		face->Y2 = this->Y + this->Width;
	}

	void SetWestInterface(Face2D* face)
	{
		this->Faces.push_back(face);
		this->WestFace = face;
		face->X1 = this->X;
		face->Y1 = this->Y;
		face->X2 = this->X;
		face->Y2 = this->Y + this->Width;
	}

	double* OuterNormalVector(Face* face)
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

	double VolumicTerm(BasisFunction<2>* phi1, BasisFunction<2>* phi2, Poisson_DG_ReferenceElement<2>* referenceElement)
	{
		return referenceElement->VolumicTerm(phi1, phi2);
		//if (this == ReferenceSquare)
			//return 
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

	function<double(Point)> EvalPhiOnFace(Face* face, BasisFunction<2>* p_phi)
	{
		IBasisFunction2D* phi = dynamic_cast<IBasisFunction2D*>(p_phi);

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


	function<double*(Point)> GradPhiOnFace(Face* face, BasisFunction<2>* p_phi)
	{
		IBasisFunction2D* phi = dynamic_cast<IBasisFunction2D*>(p_phi);

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
};

//Square* Square::ReferenceSquare = new Square(-1, -1, -1, 2);