#pragma once
#include "../../DG/Poisson_DG_Element.h"
#include "../../HHO/Poisson_HHO_Element.h"
#include "PolygonalShape.h"
using namespace std;

class Polygon : public Poisson_DG_Element<2>, public Poisson_HHO_Element<2>
{
private:
	PolygonalShape* _shape;

public:
	Polygon(int number, vector<Vertex*> vertices) :
		Element(number),
		Poisson_DG_Element<2>(number),
		Poisson_HHO_Element<2>(number)
	{
		_shape = new PolygonalShape(vertices);
	}

	inline vector<Vertex*> Vertices()
	{
		return _shape->Vertices();
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	GeometricShapeWithReferenceShape<2>* Shape() const
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

	virtual ~Polygon()
	{
		if (_shape)
			delete _shape;
	}

};