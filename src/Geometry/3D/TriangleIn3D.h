#pragma once
#include "../../Mesh/Vertex.h"
#include "../PhysicalShapeWithConstantJacobian.h"
#include "../2D/Triangle.h"
using namespace std;

class TriangleIn3D : public PhysicalShapeWithConstantJacobian<2>
{
private:
	DomPoint v1;
	DomPoint v2;
	DomPoint v3;

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
	TriangleIn3D(const DomPoint& p1, const DomPoint& p2, const DomPoint& p3)
		: v1(p1), v2(p2), v3(p3)
	{
		Init();
	}

	TriangleIn3D(const TriangleIn3D& shape) = default;

	void Init()
	{
		DimVector<3> v12 = Vect<3>(v1, v2);
		DimVector<3> v13 = Vect<3>(v1, v3);
		DimVector<3> v23 = Vect<3>(v2, v3);
		_diameter = max({ v12.norm(), v13.norm(), v23.norm() });

		_measure = 0.5 * v12.cross(v13).norm();

		_center = DomPoint((v1.X + v2.X + v3.X) / 3, (v1.Y + v2.Y + v3.Y) / 3, (v1.Z + v2.Z + v3.Z) / 3);

		_inRadius = 2 * _measure / (v12.norm() + v13.norm() + v23.norm());

		_detJacobian = _measure / RefShape()->Measure();
		
		DimMatrix<2> mapping;
		mapping <<
			v2.X - v1.X, v3.X - v1.X,
			v2.Y - v1.Y, v3.Y - v1.Y;
		if (abs(mapping.determinant()) > Utils::NumericalZero)
			_doNotUseZ = true;
		else
		{
			mapping <<
				v2.X - v1.X, v3.X - v1.X,
				v2.Z - v1.Z, v3.Z - v1.Z;
			if (abs(mapping.determinant()) > Utils::NumericalZero)
				_doNotUseY = true;
			else
			{
				mapping <<
					v2.Y - v1.Y, v3.Y - v1.Y,
					v2.Z - v1.Z, v3.Z - v1.Z;
				if (abs(mapping.determinant()) > Utils::NumericalZero)
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
	
	inline vector<DomPoint> Vertices() const override
	{
		return vector<DomPoint>{ v1, v2, v3 };
	}

	bool IsDegenerated() const override
	{
		Utils::FatalError("The function TriangleIn3D::IsDegenerated() is not implemented.");
		return true;
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
	inline DomPoint InteriorPoint() const override
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
	inline bool Contains(const DomPoint& p) const override
	{
		return Triangle::TriangleContains(v1, v2, v3, p, _measure);
	}

	void RefineWithoutCoarseOverlap(const vector<PhysicalShape<1>*>& doNotCross) override
	{
		Utils::FatalError("NOT SUPPOSED TO BE CALLED");
	}

	inline double DetJacobian() const
	{
		return _detJacobian;
	}
	inline DimMatrix<2> InverseJacobianTranspose() const
	{
		assert(false);
		return DimMatrix<2>(); // to avoid warning
	}

	DomPoint ConvertToDomain(const RefPoint& refPoint) const
	{
		double t = refPoint.X;
		double u = refPoint.Y;

		DomPoint p;
		p.X = (v2.X - v1.X) * t + (v3.X - v1.X)*u + v1.X;
		p.Y = (v2.Y - v1.Y) * t + (v3.Y - v1.Y)*u + v1.Y;
		p.Z = (v2.Z - v1.Z) * t + (v3.Z - v1.Z)*u + v1.Z;
		return p;
	}

	RefPoint ConvertToReference(const DomPoint& domainPoint) const
	{
		double x = domainPoint.X;
		double y = domainPoint.Y;
		double z = domainPoint.Z;

		DimVector<2> v;
		if (_doNotUseZ)
			v << x - v1.X, y - v1.Y;
		else if (_doNotUseY)
			v << x - v1.X, z - v1.Z;
		else if (_doNotUseX)
			v << y - v1.Y, z - v1.Z;

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
		v1.Serialize(os, 3);
		os << "--";
		v2.Serialize(os, 3);
		os << "--";
		v3.Serialize(os, 3);
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
		DomPoint lowerLeft(-1, -1, 0);
		DomPoint lowerRight(1, -1, 1);
		DomPoint upperLeft(-1, 1, 5);

		TriangleIn3D t(lowerLeft, lowerRight, upperLeft);

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
		DomPoint v1(4.5, 1.8, 0.21);
		DomPoint v2(4, 1.7, 0);
		DomPoint v3(4, 1.7, 1);
		TriangleIn3D t2(v1, v2, v3);
		t2.ConvertToReference(DomPoint(4.2, 1.7, 0.1));

	}
};