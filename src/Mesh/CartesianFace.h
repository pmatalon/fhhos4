#pragma once
#include "Vertex.h"
#include "Face.h"
#include "../Geometry/CartesianShape.h"
#include "../HHO/Diff_HHOFace.h"

template <int Dim>
class CartesianFace : public Diff_DGFace<Dim>, public Diff_HHOFace<Dim>
{
protected:
	CartesianShape<Dim, Dim - 1> _shape;
	vector<Vertex*> _vertices;
public:

	CartesianFace(BigNumber number, Vertex* origin, double width, Element<Dim>* element1, Element<Dim>* element2, CartesianShapeOrientation orientation) :
		Diff_DGFace<Dim>(number, element1, element2),
		Diff_HHOFace<Dim>(number, element1, element2),
		_shape(origin, width, orientation)
	{}

	CartesianFace(BigNumber number, Vertex* origin, double firstWidth, double secondWidth, Element<Dim>* element1, Element<Dim>* element2, CartesianShapeOrientation orientation) :
		Diff_DGFace<Dim>(number, element1, element2),
		Diff_HHOFace<Dim>(number, element1, element2),
		_shape(origin, firstWidth, secondWidth, orientation)
	{
		assert(Dim == 3);
	}

private:
	CartesianFace(BigNumber number, Element<Dim>* element1, Element<Dim>* element2, const CartesianShape<Dim, Dim - 1>& shape) :
		Face<Dim>(number, element1, element2),
		Diff_DGFace<Dim>(number, element1, element2),
		Diff_HHOFace<Dim>(number, element1, element2),
		_shape(shape)
	{}

public:

	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	const PhysicalShape<Dim-1>* Shape() const override
	{
		return &_shape;
	}
	PhysicalShape<Dim - 1>* Shape() override
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

	Face<Dim>* CreateSameGeometricFace(BigNumber number, Element<Dim>* element1)
	{
		Face<Dim>* copy = new CartesianFace<Dim>(number, element1, nullptr, _shape);
		copy->IsDomainBoundary = this->IsDomainBoundary;
		copy->BoundaryPart = this->BoundaryPart;
		return copy;
	}

	inline void ExportFaceToMatlab(FILE* file)
	{
		if (Dim == 2)
		{
			DomPoint p2 = *(_shape.Origin);
			if (_shape.Orientation == CartesianShapeOrientation::Horizontal)
				p2.X += _shape.WidthX;
			else if (_shape.Orientation == CartesianShapeOrientation::Vertical)
				p2.Y += _shape.WidthY;
			else
				assert(false);
			//             Number  x1    y1    x2    y2 IsDomainBoundary IsRemovedOnCoarserGrid
			fprintf(file, "%lu %.17g %.17g %.17g %.17g %d %d\n", this->Number, _shape.Origin->X, _shape.Origin->Y, p2.X, p2.Y, this->IsDomainBoundary, this->IsRemovedOnCoarserGrid);
		}
		else
			assert(false);
	}

	virtual ~CartesianFace()
	{}

};