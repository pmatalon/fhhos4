#pragma once
#include "../../DG/Diff_DGFace.h"
#include "../../HHO/Diff_HHOFace.h"
#include "../../Geometry/2D/Segment.h"
#include "../../Utils/Geometry.h"

class Edge : public Diff_DGFace<2>, public Diff_HHOFace<2>
{
private:
	Segment* _shape;
public:

	Edge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1, Element<2>* element2) : 
		Diff_DGFace(number, element1, element2),
		Diff_HHOFace(number, element1, element2),
		Face(number, element1, element2)
	{
		_shape = new Segment(v1, v2);
	}

	Edge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1) :
		Edge(number, v1, v2, element1, nullptr)
	{ }

	Edge(BigNumber number, Vertex* v1, Vertex* v2) :
		Edge(number, v1, v2, nullptr, nullptr)
	{ }

	inline Vertex* Vertex1() const
	{
		return _shape->Vertex1;
	}
	inline Vertex* Vertex2() const
	{
		return _shape->Vertex2;
	}
	
	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	PhysicalShape<1>* Shape() const override
	{
		return _shape;
	}

	Face<2>* CreateSameGeometricFace(BigNumber number, Element<2>* element1)
	{
		Face<2>* copy = new Edge(number, _shape->Vertex1, _shape->Vertex2, element1);
		copy->IsDomainBoundary = this->IsDomainBoundary;
		return copy;
	}

	void ExportFaceToMatlab(FILE* file)
	{
		//             Number  x1    y1    x2    y2 IsDomainBoundary IsRemovedOnCoarserGrid
		fprintf(file, "%lu %.17g %.17g %.17g %.17g %d %d\n", this->Number, _shape->Vertex1->X, _shape->Vertex1->Y, _shape->Vertex2->X, _shape->Vertex2->Y, this->IsDomainBoundary, this->IsRemovedOnCoarserGrid);
	}

	virtual ~Edge()
	{
		if (_shape)
			delete _shape;
	}
};