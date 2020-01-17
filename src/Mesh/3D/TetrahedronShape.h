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
		V1 = v1;
		V2 = v2;
		V3 = v3;
		V4 = v4;
		Init();
	}

	void Init()
	{
		double lengthEdge12 = Vect<3>(V1, V2).norm();
		double lengthEdge13 = Vect<3>(V1, V3).norm();
		double lengthEdge14 = Vect<3>(V1, V4).norm();
		double lengthEdge23 = Vect<3>(V2, V3).norm();
		double lengthEdge24 = Vect<3>(V2, V4).norm();
		double lengthEdge34 = Vect<3>(V3, V4).norm();
		_diameter = max({ lengthEdge12, lengthEdge13, lengthEdge14, lengthEdge23, lengthEdge24, lengthEdge34 });

		DimMatrix<3> m;
		m.col(0) = Vect<3>(V1, V2);
		m.col(1) = Vect<3>(V1, V3);
		m.col(2) = Vect<3>(V1, V4);
		_measure = abs(m.determinant()) / 6;

		_center = DomPoint((V1->X + V2->X + V3->X + V4->X) / 4, (V1->Y + V2->Y + V3->Y + V4->Y) / 4, (V1->Z + V2->Z + V3->Z + V4->Z) / 4);

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

	ReferenceShape<3>* RefShape() const
	{
		return &RefTetra;
	}
	
	inline vector<Vertex*> Vertices() const override
	{
		return vector<Vertex*> { V1, V2, V3, V4 };
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
	inline bool Contains(DomPoint p) const override
	{
		assert(false && "Not implemented");
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
		/*double x = domainPoint.X;
		double y = domainPoint.Y;
		double z = domainPoint.Z;*/

		DimVector<3> tuv = _inverseMapping * Vect<3>(*V1, domainPoint);

		/*double t = tuv(0);
		double u = tuv(1);
		double v = tuv(2);*/

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

		RefPoint llRef = t.ConvertToReference(A);
		assert(llRef == RefPoint(0, 0, 0));
		DomPoint llDom = t.ConvertToDomain(RefPoint(0, 0, 0));
		assert(A == llDom);

		RefPoint ulRef = t.ConvertToReference(B);
		assert(ulRef == RefPoint(1, 0, 0));

		assert(abs(t.Measure() - pow(3, 3)/6.0) < 1e-14); // Tetra's volume is 1/6 of the cube's
	}
};

ReferenceTetrahedron TetrahedronShape::RefTetra = ReferenceTetrahedron();