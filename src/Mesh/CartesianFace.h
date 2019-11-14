#pragma once
#include "Vertex.h"
#include "Face.h"
#include "CartesianShape.h"
#include "../HHO/Poisson_HHO_Face.h"

template <int Dim>
class CartesianFace : public Poisson_DG_Face<Dim>, public Poisson_HHO_Face<Dim>
{
private:
	CartesianShape<Dim, Dim - 1>* _shape;
public:

	CartesianFace(BigNumber number, Vertex* origin, double width, Element<Dim>* element1, Element<Dim>* element2, CartesianShapeOrientation orientation) :
		Poisson_DG_Face<Dim>(number, element1, element2),
		Poisson_HHO_Face<Dim>(number, element1, element2)
	{
		_shape = new CartesianShape<Dim, Dim - 1>(origin, width, orientation);
	}

	CartesianFace(BigNumber number, Vertex* origin, double firstWidth, double secondWidth, Element<Dim>* element1, Element<Dim>* element2, CartesianShapeOrientation orientation) :
		Poisson_DG_Face<Dim>(number, element1, element2),
		Poisson_HHO_Face<Dim>(number, element1, element2)
	{
		assert(Dim == 3);
		_shape = new CartesianShape<Dim, Dim - 1>(origin, firstWidth, secondWidth, orientation);
	}

private:
	CartesianFace(BigNumber number, Element<Dim>* element1, Element<Dim>* element2, CartesianShape<Dim, Dim - 1>* shape) :
		Poisson_DG_Face<Dim>(number, element1, element2),
		Poisson_HHO_Face<Dim>(number, element1, element2),
		Face<Dim>(number, element1,  element2)
	{
		_shape = shape;
	}

public:

	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	GeometricShapeWithReferenceShape<Dim-1>* Shape() const override
	{
		return _shape;
	}

	Face<Dim>* CreateSameGeometricFace(BigNumber number, Element<Dim>* element1)
	{
		Face<Dim>* copy = new CartesianFace<Dim>(number, element1, NULL, _shape);
		copy->IsDomainBoundary = this->IsDomainBoundary;
		return copy;
	}

	inline void ExportFaceToMatlab(FILE* file)
	{
		fprintf(file, "%llu %.17g %.17g %.17g %.17g %.17g %.17g %d %d\n", this->Number, _shape->Origin->X, _shape->Origin->Y, _shape->Origin->Z, _shape->WidthX, _shape->WidthY, _shape->WidthZ, _shape->Orientation, this->IsDomainBoundary);
	}

	//----------------------------------------------------------------//
	//                 Poisson_HHO_Face implementation                //
	//----------------------------------------------------------------//

	inline DenseMatrix FaceMassMatrix(FunctionalBasis<Dim-1>* basis)
	{
		return _shape->FaceMassMatrix(basis);
	}

};