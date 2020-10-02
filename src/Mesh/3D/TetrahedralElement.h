#pragma once
#include "../../DG/Diff_DGElement.h"
#include "../../HHO/Diff_HHOElement.h"
#include "../../Geometry/3D/Tetrahedron.h"
#include "TriangularFace.h"
using namespace std;

class TetrahedralElement : public Diff_DGElement<3>, public Diff_HHOElement<3>
{
private:
	Tetrahedron _shape;

	Vertex* _v1;
	Vertex* _v2;
	Vertex* _v3;
	Vertex* _v4;

public:
	TetrahedralElement(BigNumber number, Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4) :
		Element(number),
		Diff_DGElement<3>(number),
		Diff_HHOElement<3>(number),
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

	PhysicalShape<3>* Shape() override
	{
		return &_shape;
	}
	const PhysicalShape<3>* Shape() const override
	{
		return &_shape;
	}

	vector<Vertex*> Vertices() const override
	{
		return { _v1, _v2, _v3, _v4 };
	}

	DimVector<3> OuterNormalVector(Face<3>* f) const override
	{
		TriangularFace* face = dynamic_cast<TriangularFace*>(f);
		Vertex* A = face->Vertex1();
		Vertex* B = face->Vertex2();
		Vertex* C = face->Vertex3();

		DimVector<3> AB = Vect<3>(A, B);
		DimVector<3> AC = Vect<3>(A, C);

		// Condition 1: n.AB = 0 and n.AC = 0
		DimVector<3> n = AB.cross(AC);

		for (Vertex* D : this->Vertices())
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

	virtual ~TetrahedralElement()
	{}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	void UnitTests() const override
	{
		Element<3>::UnitTests();
		assert(this->Faces.size() == 4);
	}

};