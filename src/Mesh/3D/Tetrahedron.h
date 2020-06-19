#pragma once
#include "../../DG/Diff_DGElement.h"
#include "../../HHO/Diff_HHOElement.h"
#include "TetrahedronShape.h"
#include "TriangularFace.h"
using namespace std;

class Tetrahedron : public Diff_DGElement<3>, public Diff_HHOElement<3>
{
private:
	TetrahedronShape* _shape;

public:
	Tetrahedron(BigNumber number, Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4) :
		Element(number),
		Diff_DGElement<3>(number),
		Diff_HHOElement<3>(number)
	{
		_shape = new TetrahedronShape(v1, v2, v3, v4);
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

	GeometricShapeWithReferenceShape<3>* Shape() const
	{
		return _shape;
	}

	DimVector<3> OuterNormalVector(Face<3>* f) const
	{
		TriangularFace* face = dynamic_cast<TriangularFace*>(f);
		Vertex* A = face->Vertex1();
		Vertex* B = face->Vertex2();
		Vertex* C = face->Vertex3();

		DimVector<3> AB = Vect<3>(A, B);
		DimVector<3> AC = Vect<3>(A, C);

		// Condition 1: n.AB = 0 and n.AC = 0
		DimVector<3> n = AB.cross(AC);

		for (Vertex* D : this->Shape()->Vertices())
		{
			if (D != A && D != B && D != C)
			{
				// Condition 2: n.AD < 0
				DimVector<3> AD = Vect<3>(A, D);
				if (n.dot(AD) > 0)
					n = -1 * n;
			}
		}

		n = n.normalized();
		return n;
	}

	virtual ~Tetrahedron()
	{
		if (_shape)
			delete _shape;
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	void UnitTests() const override
	{
		Element<3>::UnitTests();
		assert(this->Faces.size() == 4);
	}

};