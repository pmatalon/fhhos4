#pragma once
#include "../Vertex.h"
#include "ReferenceTetrahedron.h"
#include "../GeometricShapeWithConstantJacobian.h"
using namespace std;

class TetrahedronShape : public GeometricShapeWithConstantJacobian<3>
{
private:
	double _diameter;
	double _measure;
	DomPoint _center;
	double _inRadius;

	DimMatrix<3> _inverseMapping;
	DimMatrix<3> _inverseJacobianTranspose;
	double _detJacobian;

public:
	Vertex* V1;
	Vertex* V2;
	Vertex* V3;
	Vertex* V4;

	static ReferenceTetrahedron RefTetra;

	TetrahedronShape(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4)
	{
		assert(*v1 != *v2 && *v1 != *v3 && *v1 != *v4 && *v2 != *v3 && *v2 != *v4 && *v3 != *v4);
		V1 = v1;
		V2 = v2;
		V3 = v3;
		V4 = v4;
		Init();
	}

	TetrahedronShape(const TetrahedronShape& shape) = default;

	void Init()
	{
		DimVector<3> v12 = Vect<3>(V1, V2);
		DimVector<3> v13 = Vect<3>(V1, V3);
		DimVector<3> v14 = Vect<3>(V1, V4);
		DimVector<3> v23 = Vect<3>(V2, V3);
		DimVector<3> v24 = Vect<3>(V2, V4);
		DimVector<3> v34 = Vect<3>(V3, V4);
		double lengthEdge12 = v12.norm();
		double lengthEdge13 = v13.norm();
		double lengthEdge14 = v14.norm();
		double lengthEdge23 = v23.norm();
		double lengthEdge24 = v24.norm();
		double lengthEdge34 = v34.norm();
		_diameter = max({ lengthEdge12, lengthEdge13, lengthEdge14, lengthEdge23, lengthEdge24, lengthEdge34 });

		DimMatrix<3> m;
		m.col(0) = v12;
		m.col(1) = v13;
		m.col(2) = v14;
		_measure = abs(m.determinant()) / 6;

		_center = DomPoint((V1->X + V2->X + V3->X + V4->X) / 4, (V1->Y + V2->Y + V3->Y + V4->Y) / 4, (V1->Z + V2->Z + V3->Z + V4->Z) / 4);

		_inRadius = 3 * _measure / (0.5*v12.cross(v13).norm() + 0.5*v13.cross(v14).norm() + 0.5*v14.cross(v12).norm() + 0.5*v23.cross(v24).norm());

		_detJacobian = _measure / RefTetra.Measure();


		DimMatrix<3> mapping;
		mapping <<
			V2->X - V1->X, V3->X - V1->X, V4->X - V1->X,
			V2->Y - V1->Y, V3->Y - V1->Y, V4->Y - V1->Y,
			V2->Z - V1->Z, V3->Z - V1->Z, V4->Z - V1->Z;
		_inverseMapping = mapping.inverse();

		DimMatrix<3> inverseJacobian = _inverseMapping;
		_inverseJacobianTranspose = inverseJacobian.transpose();

		//assert(abs(_detJacobian - 1 / inverseJacobian.determinant()) < 1e-14);
	}

	PhysicalShape<3>* CreateCopy() const
	{
		return new TetrahedronShape(*this);
	}

	ReferenceShape<3>* RefShape() const
	{
		return &RefTetra;
	}
	
	inline vector<Vertex*> Vertices() const override
	{
		return vector<Vertex*> { V1, V2, V3, V4 };
	}

	bool IsDegenerated() const override
	{
		assert(false && "To implement");
	}

	static ReferenceTetrahedron* InitReferenceShape()
	{
		return &RefTetra;
	}

	inline double Diameter() const override
	{
		return _diameter;
	}
	inline double Measure() const override
	{
		return _measure;
	}
	inline DomPoint Center() const override
	{
		return _center;
	}
	inline bool IsConvex() const override
	{
		return true;
	}
	inline double InRadius() const override
	{
		return _inRadius;
	}
	inline bool Contains(DomPoint p) const override
	{
		return Geometry::IsInTetrahedron(*V1, *V2, *V3, *V4, p, _measure, _diameter);
	}

	inline double DetJacobian() const
	{
		return _detJacobian;
	}
	inline DimMatrix<3> InverseJacobianTranspose() const
	{
		return _inverseJacobianTranspose;
	}

	// Mapping
	DomPoint ConvertToDomain(RefPoint refPoint) const
	{
		double t = refPoint.X;
		double u = refPoint.Y;
		double v = refPoint.Z;

		DomPoint p;
		p.X = (V2->X - V1->X)*t + (V3->X - V1->X)*u + (V4->X - V1->X)*v + V1->X;
		p.Y = (V2->Y - V1->Y)*t + (V3->Y - V1->Y)*u + (V4->Y - V1->Y)*v + V1->Y;
		p.Z = (V2->Z - V1->Z)*t + (V3->Z - V1->Z)*u + (V4->Z - V1->Z)*v + V1->Z;
		return p;
	}

	// Inverse mapping
	RefPoint ConvertToReference(DomPoint domainPoint) const
	{
		DimVector<3> tuv = _inverseMapping * Vect<3>(*V1, domainPoint);
		RefPoint p(tuv(0), tuv(1), tuv(2));
		return p;
	}

	void Serialize(ostream& os) const override
	{
		os << "Tetrahedron";
		os << " ";
		V1->Serialize(os, 3);
		os << "--";
		V2->Serialize(os, 3);
		os << "--";
		V3->Serialize(os, 3);
		os << "--";
		V4->Serialize(os, 3);
	}

	//---------------------------------------------------------------------//
	// This is f***ing useless, it should be automatic due to inheritance! //
	// But without that it doesn't compile for some reason :-(             //
	//---------------------------------------------------------------------//

	virtual double Integral(DomFunction globalFunction) const
	{
		return GeometricShapeWithConstantJacobian<3>::Integral(globalFunction);
	}
	virtual double Integral(DomFunction globalFunction, int polynomialDegree) const
	{
		return GeometricShapeWithConstantJacobian<3>::Integral(globalFunction, polynomialDegree);
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	static void Test()
	{
		int number = 0;
		Vertex A(number, 1, 1, 1);
		Vertex B(number, 4, 1, 1);
		Vertex C(number, 1, 4, 1);
		Vertex D(number, 1, 1, 4);

		TetrahedronShape t(&A, &B, &C, &D);

		t.UnitTests();

		RefPoint ARef = t.ConvertToReference(A);
		assert(ARef == RefPoint(0, 0, 0));
		DomPoint ADom = t.ConvertToDomain(RefPoint(0, 0, 0));
		assert(A == ADom);

		RefPoint BRef = t.ConvertToReference(B);
		assert(BRef == RefPoint(1, 0, 0));
		DomPoint BDom = t.ConvertToDomain(RefPoint(1, 0, 0));
		assert(B == BDom);

		RefPoint CRef = t.ConvertToReference(C);
		assert(CRef == RefPoint(0, 1, 0));
		DomPoint CDom = t.ConvertToDomain(RefPoint(0, 1, 0));
		assert(C == CDom);

		RefPoint DRef = t.ConvertToReference(D);
		assert(DRef == RefPoint(0, 0, 1));
		DomPoint DDom = t.ConvertToDomain(RefPoint(0, 0, 1));
		assert(D == DDom);

		assert(abs(t.Measure() - pow(3, 3)/6.0) < 1e-14); // Tetra's volume is 1/6 of the cube's
	}
};

ReferenceTetrahedron TetrahedronShape::RefTetra = ReferenceTetrahedron();