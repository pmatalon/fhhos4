#pragma once
#include "CartesianElement.h"
#include "IntervalFace.h"
#include "../DG/Poisson_DG_Element.h"
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

	double VolumicTerm(BasisFunction<2>* phi1, BasisFunction<2>* phi2, DiffusionPartition diffusionPartition)
	{
		double kappa = CartesianElement::DiffusionCoefficient(diffusionPartition);
		return kappa * CartesianElement::IntegralGradGrad(phi1, phi2);
	}

	double MassTerm(BasisFunction<2>* phi1, BasisFunction<2>* phi2)
	{
		return CartesianElement::MassTerm(phi1, phi2);
	}

	double SourceTerm(BasisFunction<2>* phi, SourceFunction* f)
	{
		return CartesianElement::SourceTerm(phi, f);
	}

	//-------------------------------------------------------------------//
	//                 Poisson_HHO_Element implementation                //
	//-------------------------------------------------------------------//

	Eigen::MatrixXd ComputeAndReturnCellMassMatrix(FunctionalBasis<2>* basis)
	{
		return CartesianElement::CartesianShape::ComputeAndReturnCellMassMatrix(basis);
	}

	Eigen::MatrixXd ComputeAndReturnCellReconstructMassMatrix(FunctionalBasis<2>* cellBasis, FunctionalBasis<2>* reconstructBasis)
	{
		return CartesianElement::CartesianShape::ComputeAndReturnCellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	double ComputeIntegralGradGrad(BasisFunction<2>* phi1, BasisFunction<2>* phi2)
	{
		return CartesianElement::ComputeIntegralGradGrad(phi1, phi2);
	}
	
	double St(BasisFunction<2>* reconstructPhi1, BasisFunction<2>* reconstructPhi2)
	{
		return CartesianElement::CartesianShape::IntegralGradGradReconstruct(reconstructPhi1, reconstructPhi2);
	}

	double Lt(BasisFunction<2>* phi)
	{
		return CartesianShape::Integral(phi);
	}
};