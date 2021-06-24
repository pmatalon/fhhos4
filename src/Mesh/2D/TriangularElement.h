#pragma once
#include "../../DG/Diff_DGElement.h"
#include "../../Geometry/2D/Triangle.h"
#include "Edge.h"
using namespace std;

class TriangularElement : public Diff_DGElement<2>
{
private:
	Triangle _shape;

	Vertex* _v1 = nullptr;
	Vertex* _v2 = nullptr;
	Vertex* _v3 = nullptr;

public:
	TriangularElement() {}

	TriangularElement(int number, Vertex* v1, Vertex* v2, Vertex* v3) :
		Element(number),
		Diff_DGElement<2>(number),
		_shape(*v1, *v2, *v3)
	{
		_v1 = v1;
		_v2 = v2;
		_v3 = v3;
	}

	inline Vertex* V1() { return _v1; }
	inline Vertex* V2() { return _v2; }
	inline Vertex* V3() { return _v3; }

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	PhysicalShape<2>* Shape() override
	{
		return &_shape;
	}
	const PhysicalShape<2>* Shape() const override
	{
		return &_shape;
	}

	vector<Vertex*> Vertices() const override
	{
		return { _v1, _v2, _v3 };
	}

	DimVector<2> OuterNormalVector(Face<2>* face) const override
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
		if (edge->Vertex1() == _v1)
		{
			if (edge->Vertex2() == _v2)
				C = _v3;
			else
				C = _v2;
		}
		else if (edge->Vertex1() == _v2)
		{
			if (edge->Vertex2() == _v1)
				C = _v3;
			else
				C = _v1;
		}
		else if (edge->Vertex1() == _v3)
		{
			if (edge->Vertex2() == _v1)
				C = _v2;
			else
				C = _v1;
		}
		else
			assert(false);
		
		DimVector<2> AC = Vect<2>(A, C);
		if (n.dot(AC) > 0)
			n = -1 * n;

		n = n.normalized();
		return n;
	}

	void Refine() override
	{
		_shape.RefineByConnectionOfTheMiddleEdges();
	}

	virtual ~TriangularElement()
	{}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	void UnitTests() const override
	{
		Element<2>::UnitTests();
		assert(this->Faces.size() == 3);
	}

};