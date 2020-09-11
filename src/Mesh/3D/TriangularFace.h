#pragma once
#include "../../DG/Diff_DGFace.h"
#include "../../HHO/Diff_HHOFace.h"
#include "../../Geometry/3D/TriangleIn3D.h"

class TriangularFace : public Diff_DGFace<3>, public Diff_HHOFace<3>
{
private:
	TriangleIn3D* _shape;
public:

	TriangularFace(BigNumber number, Vertex* v1, Vertex* v2, Vertex* v3, Element<3>* element1, Element<3>* element2) :
		Face(number, element1, element2),
		Diff_DGFace(number, element1, element2),
		Diff_HHOFace(number, element1, element2)
	{
		_shape = new TriangleIn3D(v1, v2, v3);
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

	PhysicalShape<2>* Shape() const override
	{
		return _shape;
	}

	Face<3>* CreateSameGeometricFace(BigNumber number, Element<3>* element1)
	{
		Face<3>* copy = new TriangularFace(number, _shape->V1, _shape->V2, _shape->V3, element1);
		copy->IsDomainBoundary = this->IsDomainBoundary;
		copy->BoundaryPart = this->BoundaryPart;
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