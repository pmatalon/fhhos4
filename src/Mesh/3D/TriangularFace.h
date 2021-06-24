#pragma once
#include "../../DG/Diff_DGFace.h"
#include "../../Geometry/3D/TriangleIn3D.h"

class TriangularFace : public Diff_DGFace<3>
{
private:
	TriangleIn3D _shape;

	Vertex* _v1;
	Vertex* _v2;
	Vertex* _v3;
public:

	TriangularFace(BigNumber number, Vertex* v1, Vertex* v2, Vertex* v3, Element<3>* element1, Element<3>* element2) :
		Face(number, element1, element2),
		Diff_DGFace(number, element1, element2),
		_shape(*v1, *v2, *v3)
	{
		_v1 = v1;
		_v2 = v2;
		_v3 = v3;
	}

	TriangularFace(BigNumber number, Vertex* v1, Vertex* v2, Vertex* v3, Element<3>* element1) :
		TriangularFace(number, v1, v2, v3, element1, nullptr)
	{}

	inline Vertex* Vertex1() const { return _v1; }
	inline Vertex* Vertex2() const { return _v2; }
	inline Vertex* Vertex3() const { return _v3; }

	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	const PhysicalShape<2>* Shape() const override
	{
		return &_shape;
	}
	PhysicalShape<2>* Shape() override
	{
		return &_shape;
	}

	vector<Vertex*> Vertices() const override
	{
		return { _v1, _v2, _v3 };
	}

	Face<3>* CreateSameGeometricFace(BigNumber number, Element<3>* element1)
	{
		Face<3>* copy = new TriangularFace(number, _v1, _v2, _v3, element1);
		copy->IsDomainBoundary = this->IsDomainBoundary;
		copy->BoundaryPart = this->BoundaryPart;
		return copy;
	}

	void ExportFaceToMatlab(FILE* file)
	{
		assert(false && "Not implemented.");
	}

	virtual ~TriangularFace()
	{}
};