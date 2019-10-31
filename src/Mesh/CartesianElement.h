#pragma once
#include "Element.h"
#include "CartesianShape.h"
#include "../HHO/Poisson_HHO_Element.h"

template <int Dim>
class CartesianElement : public Poisson_DG_Element<Dim>, public Poisson_HHO_Element<Dim>
{
protected:
	CartesianShape<Dim> _shape;

public:
	CartesianElement(BigNumber number, DomPoint* origin, double width) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number), 
		_shape(origin, width)
	{}

	CartesianElement(BigNumber number, DomPoint* origin, double widthX, double widthY) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number),
		_shape(origin, widthX, widthY)
	{}

	CartesianElement(BigNumber number, DomPoint* origin, double widthX, double widthY, double widthZ) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number),
		_shape(origin, widthX, widthY, widthZ)
	{}

	//------------------------------------------------------------------//
	//                      Element implementation                      //
	//------------------------------------------------------------------//

	const GeometricShapeWithReferenceShape<Dim>* Shape() const override
	{
		return &_shape;
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//
	
	inline double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return _shape.MassTerm(phi1, phi2);
	}

	inline double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return _shape.StiffnessTerm(phi1, phi2);
	}

	//-------------------------------------------------------------------//
	//                 Poisson_HHO_Element implementation                //
	//-------------------------------------------------------------------//

	inline DenseMatrix CellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		return _shape.CellMassMatrix(basis);
	}

	inline DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		return _shape.CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	inline double IntegralKGradGradReconstruct(Tensor<Dim>* K, BasisFunction<Dim>* reconstructPhi1, BasisFunction<Dim>* reconstructPhi2)
	{
		return _shape.IntegralKGradGradReconstruct(K, reconstructPhi1, reconstructPhi2);
	}
};