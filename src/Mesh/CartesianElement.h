#pragma once
#include "Element.h"
#include "../Geometry/CartesianShape.h"
#include "../HHO/Diff_HHOElement.h"

template <int Dim>
class CartesianElement : public Diff_DGElement<Dim>, public Diff_HHOElement<Dim>
{
protected:
	CartesianShape<Dim>* _shape;

public:
	CartesianElement(BigNumber number, DomPoint* origin, double width) :
		Element<Dim>(number),
		Diff_DGElement<Dim>(number),
		Diff_HHOElement<Dim>(number)
	{
		_shape = new CartesianShape<Dim>(origin, width);
	}

	CartesianElement(BigNumber number, DomPoint* origin, double widthX, double widthY) :
		Element<Dim>(number),
		Diff_DGElement<Dim>(number),
		Diff_HHOElement<Dim>(number)
	{
		_shape = new CartesianShape<Dim>(origin, widthX, widthY);
	}

	CartesianElement(BigNumber number, DomPoint* origin, double widthX, double widthY, double widthZ) :
		Element<Dim>(number),
		Diff_DGElement<Dim>(number),
		Diff_HHOElement<Dim>(number)
	{
		_shape = new CartesianShape<Dim>(origin, widthX, widthY, widthZ);
	}

	//------------------------------------------------------------------//
	//                      Element implementation                      //
	//------------------------------------------------------------------//

	PhysicalShape<Dim>* Shape() const override
	{
		return _shape;
	}

	virtual ~CartesianElement()
	{
		if (_shape)
			delete _shape;
	}
};