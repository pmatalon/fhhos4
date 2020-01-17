#pragma once
#include "../../DG/Poisson_DG_Face.h"
#include "../../HHO/Poisson_HHO_Face.h"
#include "TriangleIn3DShape.h"

class TriangularFace : public Poisson_DG_Face<3>, public Poisson_HHO_Face<3>
{
private:
	TriangleIn3DShape* _shape;
public:

	TriangularFace(BigNumber number, Vertex* v1, Vertex* v2, Vertex* v3, Element<3>* element1, Element<3>* element2) :
		Poisson_DG_Face(number, element1, element2),
		Poisson_HHO_Face(number, element1, element2),
		Face(number, element1, element2)
	{
		_shape = new TriangleIn3DShape(v1, v2, v3);
	}

	TriangularFace(BigNumber number, Vertex* v1, Vertex* v2, Vertex* v3, Element<3>* element1) :
		TriangularFace(number, v1, v2, v3, element1, nullptr)
	{ }

	inline Vertex* Vertex1() const
	{
		return _shape->V1;
	}
	inline Vertex* Vertex2() const
	{
		return _shape->V2;
	}
	inline Vertex* Vertex3() const
	{
		return _shape->V3;
	}

	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	GeometricShapeWithReferenceShape<2>* Shape() const override
	{
		return _shape;
	}

	Face<3>* CreateSameGeometricFace(BigNumber number, Element<3>* element1)
	{
		Face<3>* copy = new TriangularFace(number, _shape->V1, _shape->V2, _shape->V3, element1);
		copy->IsDomainBoundary = this->IsDomainBoundary;
		return copy;
	}

	void ExportFaceToMatlab(FILE* file)
	{
		assert(false && "Not implemented.");
	}

	virtual ~TriangularFace()
	{
		if (_shape)
			delete _shape;
	}
};