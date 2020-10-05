#pragma once
#include "../../DG/Diff_DGFace.h"
#include "../../HHO/Diff_HHOFace.h"
#include "../../Geometry/2D/Segment.h"

class Edge : public Diff_DGFace<2>, public Diff_HHOFace<2>
{
private:
	Segment _shape;
	Vertex* _v1 = nullptr;
	Vertex* _v2 = nullptr;
public:
	Edge() {}

	Edge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1, Element<2>* element2) :
		Face(number, element1, element2),
		Diff_DGFace(number, element1, element2),
		Diff_HHOFace(number, element1, element2),
		_shape(*v1, *v2)
	{
		_v1 = v1;
		_v2 = v2;
	}

	Edge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1) :
		Edge(number, v1, v2, element1, nullptr)
	{}

	Edge(BigNumber number, Vertex* v1, Vertex* v2) :
		Edge(number, v1, v2, nullptr, nullptr)
	{}

	inline Vertex* Vertex1() const
	{
		return _v1;
	}
	inline Vertex* Vertex2() const
	{
		return _v2;
	}
	
	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	const PhysicalShape<1>* Shape() const override
	{
		return &_shape;
	}
	PhysicalShape<1>* Shape() override
	{
		return &_shape;
	}

	vector<Vertex*> Vertices() const override
	{
		return { _v1, _v2 };
	}

	Face<2>* CreateSameGeometricFace(BigNumber number, Element<2>* element1)
	{
		Face<2>* copy = new Edge(number, _v1, _v2, element1);
		copy->IsDomainBoundary = this->IsDomainBoundary;
		copy->BoundaryPart = this->BoundaryPart;
		return copy;
	}

	void ExportFaceToMatlab(FILE* file)
	{
		//             Number  x1    y1    x2    y2 IsDomainBoundary IsRemovedOnCoarserGrid
		fprintf(file, "%lu %.17g %.17g %.17g %.17g %d %d\n", this->Number, _v1->X, _v1->Y, _v2->X, _v2->Y, this->IsDomainBoundary, this->IsRemovedOnCoarserGrid);
	}

	virtual ~Edge()
	{}
};