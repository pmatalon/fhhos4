#pragma once
#include "Element.h"
#include "CartesianShape.h"
#include "../HHO/Poisson_HHO_Element.h"

template <int Dim>
class CartesianElement : public Poisson_DG_Element<Dim>, public Poisson_HHO_Element<Dim>
{
protected:
	CartesianShape<Dim>* _shape;

public:
	CartesianElement(BigNumber number, DomPoint* origin, double width) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number)
	{
		_shape = new CartesianShape<Dim>(origin, width);
	}

	CartesianElement(BigNumber number, DomPoint* origin, double widthX, double widthY) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number)
	{
		_shape = new CartesianShape<Dim>(origin, widthX, widthY);
	}

	CartesianElement(BigNumber number, DomPoint* origin, double widthX, double widthY, double widthZ) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number)
	{
		_shape = new CartesianShape<Dim>(origin, widthX, widthY, widthZ);
	}

	//------------------------------------------------------------------//
	//                      Element implementation                      //
	//------------------------------------------------------------------//

	GeometricShapeWithReferenceShape<Dim>* Shape() const override
	{
		return _shape;
	}

	virtual ~CartesianElement()
	{
		if (_shape)
			delete _shape;
	}
};