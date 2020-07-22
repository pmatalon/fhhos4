#pragma once
#include "../../DG/Diff_DGElement.h"
#include "../../HHO/Diff_HHOElement.h"
#include "TriangleShape.h"
#include "Edge.h"
using namespace std;

class Triangle : public Diff_DGElement<2>, public Diff_HHOElement<2>
{
private:
	TriangleShape* _shape;

public:
	Triangle(int number, Vertex* v1, Vertex* v2, Vertex* v3) :
		Element(number),
		Diff_DGElement<2>(number),
		Diff_HHOElement<2>(number)
	{
		_shape = new TriangleShape(v1, v2, v3);
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
		Vertex* C = nullptr;
		if (edge->Vertex1() == _shape->V1)
		{
			if (edge->Vertex2() == _shape->V2)
				C = _shape->V3;
			else
				C = _shape->V2;
		}
		else if (edge->Vertex1() == _shape->V2)
		{
			if (edge->Vertex2() == _shape->V1)
				C = _shape->V3;
			else
				C = _shape->V1;
		}
		else if (edge->Vertex1() == _shape->V3)
		{
			if (edge->Vertex2() == _shape->V1)
				C = _shape->V2;
			else
				C = _shape->V1;
		}
		else
			assert(false);
		
		DimVector<2> AC = Vect<2>(A, C);
		if (n.dot(AC) > 0)
			n = -1 * n;

		n = n.normalized();
		return n;
	}

	virtual ~Triangle()
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
		assert(this->Faces.size() == 3);
	}

};