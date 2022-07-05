#pragma once
#include "Element.h"
#include "../Geometry/CartesianShape.h"
#include "../Discretizations/DG/Diff_DGElement.h"

template <int Dim>
class CartesianElement : public Diff_DGElement<Dim>
{
protected:
	CartesianShape<Dim> _shape;
	vector<Vertex*> _vertices;
public:
	CartesianElement() {}

	CartesianElement(BigNumber number, DomPoint* origin, double width) :
		Element<Dim>(number),
		Diff_DGElement<Dim>(number),
		_shape(origin, width)
	{}

	CartesianElement(BigNumber number, DomPoint* origin, double widthX, double widthY) :
		Element<Dim>(number),
		Diff_DGElement<Dim>(number),
		_shape(origin, widthX, widthY)
	{}

	CartesianElement(BigNumber number, DomPoint* origin, double widthX, double widthY, double widthZ) :
		Element<Dim>(number),
		Diff_DGElement<Dim>(number),
		_shape(origin, widthX, widthY, widthZ)
	{}

	//------------------------------------------------------------------//
	//                      Element implementation                      //
	//------------------------------------------------------------------//

	PhysicalShape<Dim>* Shape() override
	{
		return &_shape;
	}
	const PhysicalShape<Dim>* Shape() const override
	{
		return &_shape;
	}

	inline void SetVertices(const vector<Vertex*>& vertices)
	{
		_vertices = vertices;
		_shape.SetVertices(Vertex::ToDomPoints(vertices));
	}

	vector<Vertex*> Vertices() const override
	{
		return _vertices;
	}

	void Refine(int nRefinements) override
	{
		Utils::FatalError("Method CartesianElement::Refine() is not implemented!");
	}

	virtual ~CartesianElement()
	{}
};