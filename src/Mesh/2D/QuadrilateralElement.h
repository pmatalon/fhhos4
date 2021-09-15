#pragma once
#include "../../DG/Diff_DGElement.h"
#include "../../Geometry/2D/Quadrilateral.h"
#include "Edge.h"
using namespace std;

class QuadrilateralElement : public Diff_DGElement<2>
{
private:
	Quadrilateral _shape;

	Vertex* _v1 = nullptr;
	Vertex* _v2 = nullptr;
	Vertex* _v3 = nullptr;
	Vertex* _v4 = nullptr;

public:
	QuadrilateralElement() {}

	QuadrilateralElement(int number, Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4) :
		Element(number),
		Diff_DGElement<2>(number),
		_shape(*v1, *v2, *v3, *v4)
	{
		_v1 = v1;
		_v2 = v2;
		_v3 = v3;
		_v4 = v4;
	}

	inline Vertex* V1() { return _v1; }
	inline Vertex* V2() { return _v2; }
	inline Vertex* V3() { return _v3; }
	inline Vertex* V4() { return _v4; }

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
		return { _v1, _v2, _v3, _v4 };
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

	void Refine(int nRefinements) override
	{
		Utils::FatalError("Method QuadrilateralElement::Refine() is not implemented!");
	}
	
	virtual ~QuadrilateralElement()
	{}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	void UnitTests() const override
	{
		Element<2>::UnitTests();
		assert(this->Faces.size() == 4);
	}

};