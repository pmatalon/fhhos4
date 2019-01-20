#pragma once
#include "Element.h"
#include "Poisson_DG_Element.h"
#include "Poisson_DG_ReferenceElement.h"

class Interval : public Element, public Poisson_DG_Element
{
public:
	double A;
	double B;

	Face* Left;
	Face* Right;

	Interval(BigNumber number, double a, double b, Face* left, Face* right) : Element(number)
	{
		this->A = a;
		this->B = b;
		this->Faces.push_back(left);
		this->Faces.push_back(right);
		this->Left = left;
		this->Right = right;
	}

	StandardElementCode StdElementCode()
	{
		return StandardElementCode::Interval;
	}

	double* OuterNormalVector(Face* interface)
	{
		if (interface == this->Left)
			return new double[1]{ -1 };
		else if(interface == this->Right)
			return new double[1]{ 1 };
		return NULL;
	}
	
	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//

	double VolumicTerm(BasisFunction* phi1, BasisFunction* phi2, Poisson_DG_ReferenceElement* referenceElement)
	{
		double h = this->B - this->A;
		double factor = 2 / h;
		return factor * referenceElement->VolumicTerm(phi1, phi2);
	}

	function<double(Point)> EvalPhiOnFace(Face* face, BasisFunction* p_phi)
	{
		IBasisFunction1D* phi = dynamic_cast<IBasisFunction1D*>(p_phi);

		double tFixed = face == this->Left ? -1 : 1;
		function<double(Point)> evalOnFace = [phi, tFixed](Point point0D) {
			return phi->Eval(tFixed);
		};
		return evalOnFace;
	}


	function<double*(Point)> GradPhiOnFace(Face* face, BasisFunction* p_phi)
	{
		IBasisFunction1D* phi = dynamic_cast<IBasisFunction1D*>(p_phi);

		double tFixed = face == this->Left ? -1 : 1;
		function<double*(Point)> gradOnFace = [phi, tFixed](Point point0D) {
			return phi->Grad(tFixed);
		};
		return gradOnFace;
	}
};