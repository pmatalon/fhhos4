#pragma once
#include "Face.h"
#include "CartesianShape.h"
#include "../HHO/Poisson_HHO_Face.h"

template <int Dim>
class CartesianFace : public Poisson_DG_Face<Dim>, public Poisson_HHO_Face<Dim>, public CartesianShape<Dim, Dim-1>
{
public:

	CartesianFace(BigNumber number, DomPoint origin, double width, Element<Dim>* element1, Element<Dim>* element2, CartesianShapeOrientation orientation) :
		Poisson_DG_Face<Dim>(number, element1, element2),
		Poisson_HHO_Face<Dim>(number, element1, element2), 
		CartesianShape<Dim, Dim - 1>(origin, width, orientation)
	{}

	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	double GetDiameter()
	{
		return CartesianShape<Dim, Dim - 1>::Width;
	}

	double Measure()
	{
		return CartesianShape<Dim, Dim - 1>::Measure();
	}

	RefPoint ConvertToReference(DomPoint domainPoint)
	{
		return CartesianShape<Dim, Dim - 1>::ConvertToReference(domainPoint);
	}
	DomPoint ConvertToDomain(RefPoint referenceElementPoint)
	{
		return CartesianShape<Dim, Dim - 1>::ConvertToDomain(referenceElementPoint);
	}

	double ComputeIntegral(function<double(RefPoint)> func, int numberOfDerivatives)
	{
		return CartesianShape<Dim, Dim - 1>::ComputeIntegral(func, numberOfDerivatives);
	}
	double ComputeIntegral(function<double(RefPoint)> func, int numberOfDerivatives, int polynomialDegree)
	{
		return CartesianShape<Dim, Dim - 1>::ComputeIntegral(func, numberOfDerivatives, polynomialDegree);
	}

	//----------------------------------------------------------------//
	//                 Poisson_HHO_Face implementation                //
	//----------------------------------------------------------------//

	Eigen::MatrixXd FaceMassMatrix(FunctionalBasis<Dim-1>* basis)
	{
		return CartesianShape<Dim, Dim - 1>::FaceMassMatrix(basis);
	}

};