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

	inline double GetDiameter() override
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
	
	inline double IntegralGlobalFunction(DomFunction func) const override
	{
		return CartesianShape<Dim>::IntegralGlobalFunction(func);
	}

	inline double Integral(BasisFunction<Dim>* phi) const
	{
		return CartesianShape<Dim>::Integral(phi);
	}

	inline double ComputeIntegral(RefFunction func) const
	{
		return CartesianShape<Dim>::ComputeIntegral(func);
	}

	inline double ComputeIntegral(RefFunction func, int polynomialDegree) const
	{
		return CartesianShape<Dim>::ComputeIntegral(func, polynomialDegree);
	}

	inline double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return CartesianShape<Dim>::ComputeIntegralGradGrad(phi1, phi2);
	}

	inline double ComputeIntegralKGradGrad(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return CartesianShape<Dim>::ComputeIntegralKGradGrad(K, phi1, phi2);
	}

	inline DomPoint ConvertToDomain(RefPoint referenceElementPoint) const
	{
		return CartesianShape<Dim>::ConvertToDomain(referenceElementPoint);
	}

	inline RefPoint ConvertToReference(DomPoint domainPoint) const
	{
		return CartesianShape<Dim>::ConvertToReference(domainPoint);
	}

	inline DimMatrix<Dim> InverseJacobian() const
	{
		return CartesianShape<Dim>::InverseJacobian();
	}

	inline vector<RefPoint> GetNodalPoints(FunctionalBasis<Dim>* basis) const
	{
		return CartesianShape<Dim>::GetNodalPoints(basis);
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//
	
	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return CartesianShape<Dim>::MassTerm(phi1, phi2);
	}

	double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return CartesianShape<Dim>::StiffnessTerm(phi1, phi2);
	}

	//-------------------------------------------------------------------//
	//                 Poisson_HHO_Element implementation                //
	//-------------------------------------------------------------------//

	DenseMatrix CellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		return CartesianShape<Dim>::CellMassMatrix(basis);
	}

	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		return CartesianShape<Dim>::CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	double IntegralKGradGradReconstruct(Tensor<Dim>* K, BasisFunction<Dim>* reconstructPhi1, BasisFunction<Dim>* reconstructPhi2)
	{
		return CartesianShape<Dim>::IntegralKGradGradReconstruct(K, reconstructPhi1, reconstructPhi2);
	}
};