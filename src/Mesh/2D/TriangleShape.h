#pragma once
#include "../Vertex.h"
#include "../ReferenceTriangle.h"
#include "../GeometricShapeWithConstantJacobian.h"
using namespace std;

class TriangleShape : public GeometricShapeWithConstantJacobian<2>
{
private:
	double _diameter;
	double _measure;
	DomPoint _center;
	double _inRadius;

	DimMatrix<2> _inverseJacobianTranspose;
	double _detJacobian;

public:
	Vertex* V1;
	Vertex* V2;
	Vertex* V3;

	static ReferenceTriangle RefTriangle;

	TriangleShape(Vertex* v1, Vertex* v2, Vertex* v3) 
	{
		V1 = v1;
		V2 = v2;
		V3 = v3;
		Init();
	}

	void Init()
	{
		double lengthEdge12 = sqrt(pow(V2->X - V1->X, 2) + pow(V2->Y - V1->Y, 2));
		double lengthEdge23 = sqrt(pow(V3->X - V2->X, 2) + pow(V3->Y - V2->Y, 2));
		double lengthEdge13 = sqrt(pow(V3->X - V1->X, 2) + pow(V3->Y - V1->Y, 2));
		_diameter = max(lengthEdge12, max(lengthEdge23, lengthEdge13));

		_measure = 0.5 * abs(V1->X * (V2->Y - V3->Y) + V2->X * (V3->Y - V1->Y) + V3->X * (V1->Y - V2->Y));

		_inRadius = 2 * _measure / (lengthEdge12 + lengthEdge23 + lengthEdge13);

		_center = DomPoint((V1->X + V2->X + V3->X) / 3, (V1->Y + V2->Y + V3->Y) / 3);

		_detJacobian = _measure / RefTriangle.Measure();

		DimMatrix<2> inverseJacobian;
		inverseJacobian(0, 0) = (V3->Y - V1->Y) / ((V3->Y - V1->Y)*(V2->X - V1->X) - (V3->X - V1->X)*(V2->Y - V1->Y));
		inverseJacobian(0, 1) = -(V3->X - V1->X) / ((V3->Y - V1->Y)*(V2->X - V1->X) - (V3->X - V1->X)*(V2->Y - V1->Y));
		inverseJacobian(1, 0) = (V2->Y - V1->Y) / ((V2->Y - V1->Y)*(V3->X - V1->X) - (V2->X - V1->X)*(V3->Y - V1->Y));
		inverseJacobian(1, 1) = -(V2->X - V1->X) / ((V2->Y - V1->Y)*(V3->X - V1->X) - (V2->X - V1->X)*(V3->Y - V1->Y));
		_inverseJacobianTranspose = inverseJacobian.transpose();
	}

	ReferenceShape<2>* RefShape() const
	{
		return &RefTriangle;
	}
	
	inline vector<Vertex*> Vertices() const override
	{
		return vector<Vertex*> { V1, V2, V3 };
	}

	static ReferenceTriangle* InitReferenceShape()
	{
		return &RefTriangle;
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
	inline double InRadius() const override
	{
		return _inRadius;
	}
	inline bool Contains(DomPoint p) const override
	{
		assert(false && "Not implemented");
	}

	inline double DetJacobian() const
	{
		return _detJacobian;
	}
	inline DimMatrix<2> InverseJacobianTranspose() const
	{
		return _inverseJacobianTranspose;
	}

	DomPoint ConvertToDomain(RefPoint refPoint) const
	{
		double t = refPoint.X;
		double u = refPoint.Y;

		DomPoint p;
		p.X = (V2->X - V1->X) * t + (V3->X - V1->X)*u + V1->X;
		p.Y = (V2->Y - V1->Y) * t + (V3->Y - V1->Y)*u + V1->Y;
		return p;
	}

	RefPoint ConvertToReference(DomPoint domainPoint) const
	{
		double x = domainPoint.X;
		double y = domainPoint.Y;

		double t = ((V3->Y - V1->Y)*(x - V1->X) - (V3->X - V1->X)*(y - V1->Y)) / ((V3->Y - V1->Y)*(V2->X - V1->X) - (V3->X - V1->X)*(V2->Y - V1->Y));
		double u = ((V2->Y - V1->Y)*(x - V1->X) - (V2->X - V1->X)*(y - V1->Y)) / ((V2->Y - V1->Y)*(V3->X - V1->X) - (V2->X - V1->X)*(V3->Y - V1->Y));
		RefPoint p(t, u);
		return p;
	}

	void Serialize(ostream& os) const override
	{
		os << "Triangle";
		os << " ";
		V1->Serialize(os, 2);
		os << "--";
		V2->Serialize(os, 2);
		os << "--";
		V3->Serialize(os, 2);
	}

	//---------------------------------------------------------------------//
	// This is f***ing useless, it should be automatic due to inheritance! //
	// But without that it doesn't compile for some reason :-(             //
	//---------------------------------------------------------------------//

	virtual double Integral(DomFunction globalFunction) const
	{
		return GeometricShapeWithConstantJacobian<2>::Integral(globalFunction);
	}
	virtual double Integral(DomFunction globalFunction, int polynomialDegree) const
	{
		return GeometricShapeWithConstantJacobian<2>::Integral(globalFunction, polynomialDegree);
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	static void Test()
	{
		int number = 0;
		Vertex lowerLeft(number, -1, -1);
		Vertex lowerRight(number, 1, -1);
		Vertex upperLeft(number, -1, 1);

		TriangleShape t(&lowerLeft, &lowerRight, &upperLeft);

		t.UnitTests();

		RefPoint llRef = t.ConvertToReference(lowerLeft);
		assert(llRef == RefPoint(0, 0));
		DomPoint llDom = t.ConvertToDomain(RefPoint(0, 0));
		assert(lowerLeft == llDom);

		RefPoint ulRef = t.ConvertToReference(upperLeft);
		assert(ulRef == RefPoint(0, 1));
	}
};

ReferenceTriangle TriangleShape::RefTriangle = ReferenceTriangle();