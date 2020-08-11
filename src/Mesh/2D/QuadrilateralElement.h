#pragma once
#include "../../DG/Diff_DGElement.h"
#include "../../HHO/Diff_HHOElement.h"
#include "../../Geometry/2D/QuadrilateralShape.h"
#include "Edge.h"
using namespace std;

class QuadrilateralElement : public Diff_DGElement<2>, public Diff_HHOElement<2>
{
private:
	QuadrilateralShape* _shape;

public:
	QuadrilateralElement(int number, Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4) :
		Element(number),
		Diff_DGElement<2>(number),
		Diff_HHOElement<2>(number)
	{
		_shape = new QuadrilateralShape(v1, v2, v3, v4);
	}

	inline Vertex* V1()
	{
		return _shape->V1;
	}
	inline Vertex* V2()
	{
		return _shape->V2;
	}
	inline Vertex* V3()
	{
		return _shape->V3;
	}
	inline Vertex* V4()
	{
		return _shape->V4;
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	PhysicalShape<2>* Shape() const
	{
		return _shape;
	}

	DimVector<2> OuterNormalVector(Face<2>* face) const
	{
		DimVector<2> n;
		Edge* edge = dynamic_cast<Edge*>(face);
		Vertex* A = edge->Vertex1();
		Vertex* B = edge->Vertex2();

		// Condition 1: n.AB = 0
		// =>  n = (-AB.Y, AB.X)
		n << A->Y - B->Y, B->X - A->X;

		// Condition 2: n.AC < 0
		DimVector<2> AC = this->Center() - *A;
		if (n.dot(AC) > 0)
			n = -1 * n;

		n = n.normalized();
		return n;
	}
	
	virtual ~QuadrilateralElement()
	{
		if (_shape)
			delete _shape;
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	void UnitTests() const override
	{
		Element<2>::UnitTests();
		assert(this->Faces.size() == 4);
	}

};