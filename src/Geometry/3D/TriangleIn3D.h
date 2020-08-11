#pragma once
#include "../../Mesh/Vertex.h"
#include "../PhysicalShapeWithConstantJacobian.h"
#include "../2D/Triangle.h"
using namespace std;

class TriangleIn3D : public PhysicalShapeWithConstantJacobian<2>
{
private:
	double _diameter;
	double _measure;
	DomPoint _center;
	double _inRadius;

	DimMatrix<2> _inverseMapping;
	bool _doNotUseX = false;
	bool _doNotUseY = false;
	bool _doNotUseZ = false;

	double _detJacobian;

public:
	Vertex* V1;
	Vertex* V2;
	Vertex* V3;

	TriangleIn3D(Vertex* v1, Vertex* v2, Vertex* v3)
	{
		V1 = v1;
		V2 = v2;
		V3 = v3;
		Init();
	}

	TriangleIn3D(const TriangleIn3D& shape) = default;

	void Init()
	{
		DimVector<3> v12 = Vect<3>(V1, V2);
		DimVector<3> v13 = Vect<3>(V1, V3);
		DimVector<3> v23 = Vect<3>(V2, V3);
		_diameter = max({ v12.norm(), v13.norm(), v23.norm() });

		_measure = 0.5 * v12.cross(v13).norm();

		_center = DomPoint((V1->X + V2->X + V3->X) / 3, (V1->Y + V2->Y + V3->Y) / 3, (V1->Z + V2->Z + V3->Z) / 3);

		_inRadius = 2 * _measure / (v12.norm() + v13.norm() + v23.norm());

		_detJacobian = _measure / RefShape()->Measure();
		
		DimMatrix<2> mapping;
		mapping <<
			V2->X - V1->X, V3->X - V1->X,
			V2->Y - V1->Y, V3->Y - V1->Y;
		if (abs(mapping.determinant()) > 1e-12)
			_doNotUseZ = true;
		else
		{
			mapping <<
				V2->X - V1->X, V3->X - V1->X,
				V2->Z - V1->Z, V3->Z - V1->Z;
			if (abs(mapping.determinant()) > 1e-12)
				_doNotUseY = true;
			else
			{
				mapping <<
					V2->Y - V1->Y, V3->Y - V1->Y,
					V2->Z - V1->Z, V3->Z - V1->Z;
				if (abs(mapping.determinant()) > 1e-12)
					_doNotUseX = true;
				else
					assert(false);
			}
		}
		_inverseMapping = mapping.inverse();
	}

	PhysicalShape<2>* CreateCopy() const
	{
		return new TriangleIn3D(*this);
	}

	ReferenceShape<2>* RefShape() const
	{
		return &Triangle::RefTriangle;
	}
	
	inline vector<Vertex*> Vertices() const override
	{
		return vector<Vertex*> { V1, V2, V3 };
	}

	bool IsDegenerated() const override
	{
		assert(false && "To implement");
	}

	static ReferenceTriangle* InitReferenceShape()
	{
		return &Triangle::RefTriangle;
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
		return Geometry::IsInTriangle(*V1, *V2, *V3, p, _measure, _diameter);
	}

	inline double DetJacobian() const
	{
		return _detJacobian;
	}
	inline DimMatrix<2> InverseJacobianTranspose() const
	{
		assert(false);
	}

	DomPoint ConvertToDomain(RefPoint refPoint) const
	{
		double t = refPoint.X;
		double u = refPoint.Y;

		DomPoint p;
		p.X = (V2->X - V1->X) * t + (V3->X - V1->X)*u + V1->X;
		p.Y = (V2->Y - V1->Y) * t + (V3->Y - V1->Y)*u + V1->Y;
		p.Z = (V2->Z - V1->Z) * t + (V3->Z - V1->Z)*u + V1->Z;
		return p;
	}

	RefPoint ConvertToReference(DomPoint domainPoint) const
	{
		double x = domainPoint.X;
		double y = domainPoint.Y;
		double z = domainPoint.Z;

		DimVector<2> v;
		if (_doNotUseZ)
			v << x - V1->X, y - V1->Y;
		else if (_doNotUseY)
			v << x - V1->X, z - V1->Z;
		else if (_doNotUseX)
			v << y - V1->Y, z - V1->Z;

		DimVector<2> tu = _inverseMapping * v;
		double t = tu(0);
		double u = tu(1);

		//if (abs(t) > 1.1 || abs(u) > 1.1 || t != t || u != u)
			//assert(false);

		RefPoint p(t, u);
		return p;
	}

	void Serialize(ostream& os) const override
	{
		os << "Triangle";
		os << " ";
		V1->Serialize(os, 3);
		os << "--";
		V2->Serialize(os, 3);
		os << "--";
		V3->Serialize(os, 3);
	}

	//---------------------------------------------------------------------//
	// This is f***ing useless, it should be automatic due to inheritance! //
	// But without that it doesn't compile for some reason :-(             //
	//---------------------------------------------------------------------//

	virtual double Integral(DomFunction globalFunction) const
	{
		return PhysicalShapeWithConstantJacobian<2>::Integral(globalFunction);
	}
	virtual double Integral(DomFunction globalFunction, int polynomialDegree) const
	{
		return PhysicalShapeWithConstantJacobian<2>::Integral(globalFunction, polynomialDegree);
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	static void Test()
	{
		int number = 0;
		Vertex lowerLeft(number, -1, -1, 0);
		Vertex lowerRight(number, 1, -1, 1);
		Vertex upperLeft(number, -1, 1, 5);

		TriangleIn3D t(&lowerLeft, &lowerRight, &upperLeft);

		t.UnitTests();

		RefPoint llRef = t.ConvertToReference(lowerLeft);
		assert(llRef == RefPoint(0, 0));
		DomPoint llDom = t.ConvertToDomain(RefPoint(0, 0, 0));
		assert(lowerLeft == llDom);

		RefPoint lrRef = t.ConvertToReference(lowerRight);
		assert(lrRef == RefPoint(1, 0));
		DomPoint lrDom = t.ConvertToDomain(RefPoint(1, 0, 0));
		assert(lowerRight == lrDom);

		RefPoint ulRef = t.ConvertToReference(upperLeft);
		assert(ulRef == RefPoint(0, 1));
		DomPoint ulDom = t.ConvertToDomain(RefPoint(0, 1, 0));
		assert(upperLeft == ulDom);

		//-----------------------------------------------------
		Vertex v1(number, 4.5, 1.8, 0.21);
		Vertex v2(number, 4, 1.7, 0);
		Vertex v3(number, 4, 1.7, 1);
		TriangleIn3D t2(&v1, &v2, &v3);
		t2.ConvertToReference(DomPoint(4.2, 1.7, 0.1));

	}
};