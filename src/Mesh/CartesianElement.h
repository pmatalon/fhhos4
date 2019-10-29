#pragma once
#include "Element.h"
#include "CartesianShape.h"
#include "../HHO/Poisson_HHO_Element.h"

template <int Dim>
class CartesianElement : public Poisson_DG_Element<Dim>, public Poisson_HHO_Element<Dim>, public CartesianShape<Dim>
{
public:
	CartesianElement(BigNumber number, DomPoint* origin, double width) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number), 
		CartesianShape<Dim>(origin, width)
	{}

	CartesianElement(BigNumber number, DomPoint* origin, double widthX, double widthY) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number),
		CartesianShape<Dim>(origin, widthX, widthY)
	{}

	CartesianElement(BigNumber number, DomPoint* origin, double widthX, double widthY, double widthZ) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number),
		CartesianShape<Dim>(origin, widthX, widthY, widthZ)
	{}

	void Serialize(ostream& os) const override
	{
		Element<Dim>::Serialize(os);
		os << ", ";
		CartesianShape<Dim>::Serialize(os);
	}

	//------------------------------------------------------------------//
	//                      Element implementation                      //
	//------------------------------------------------------------------//

	// Geometric information
	inline double Diameter() override
	{
		return max({ CartesianShape<Dim>::WidthX, CartesianShape<Dim>::WidthY, CartesianShape<Dim>::WidthZ });
	}
	inline double Measure()
	{
		return CartesianShape<Dim>::Measure;
	}
	inline DomPoint Center()
	{
		return CartesianShape<Dim>::Center;
	}

	// Transformation to reference element
	inline DomPoint ConvertToDomain(RefPoint referenceElementPoint) const
	{
		return CartesianShape<Dim>::ConvertToDomain(referenceElementPoint);
	}
	inline RefPoint ConvertToReference(DomPoint domainPoint) const
	{
		return CartesianShape<Dim>::ConvertToReference(domainPoint);
	}
	inline DimMatrix<Dim> InverseJacobianTranspose() const
	{
		return CartesianShape<Dim>::InverseJacobianTranspose();
	}
	
	// Integral
	inline double Integral(DomFunction func) const override
	{
		return CartesianShape<Dim>::Integral(func);
	}
	inline double Integral(BasisFunction<Dim>* phi) const
	{
		return CartesianShape<Dim>::Integral(phi);
	}
	inline double Integral(RefFunction func) const
	{
		return CartesianShape<Dim>::Integral(func);
	}
	inline double Integral(RefFunction func, int polynomialDegree) const
	{
		return CartesianShape<Dim>::Integral(func, polynomialDegree);
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//
	
	inline double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return CartesianShape<Dim>::MassTerm(phi1, phi2);
	}

	inline double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return CartesianShape<Dim>::StiffnessTerm(phi1, phi2);
	}

	//-------------------------------------------------------------------//
	//                 Poisson_HHO_Element implementation                //
	//-------------------------------------------------------------------//

	inline DenseMatrix CellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		return CartesianShape<Dim>::CellMassMatrix(basis);
	}

	inline DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		return CartesianShape<Dim>::CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	inline double IntegralKGradGradReconstruct(Tensor<Dim>* K, BasisFunction<Dim>* reconstructPhi1, BasisFunction<Dim>* reconstructPhi2)
	{
		return CartesianShape<Dim>::IntegralKGradGradReconstruct(K, reconstructPhi1, reconstructPhi2);
	}

	inline double ComputeIntegralKGradGrad(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return CartesianShape<Dim>::ComputeIntegralKGradGrad(K, phi1, phi2);
	}
};