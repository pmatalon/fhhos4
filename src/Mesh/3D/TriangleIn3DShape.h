#pragma once
#include "../Vertex.h"
#include "../ReferenceTriangle.h"
#include "../GeometricShapeWithConstantJacobian.h"
#include "../2D/TriangleShape.h"
using namespace std;

class TriangleIn3DShape : public GeometricShapeWithConstantJacobian<2>
{
private:
	double _diameter;
	double _measure;
	DomPoint _center;

	//DimMatrix<2> _inverseJacobianTranspose;
	double _detJacobian;

public:
	Vertex* V1;
	Vertex* V2;
	Vertex* V3;

	//static ReferenceTriangle RefTriangle;

	TriangleIn3DShape(Vertex* v1, Vertex* v2, Vertex* v3)
	{
		V1 = v1;
		V2 = v2;
		V3 = v3;
		Init();
	}

	void Init()
	{
		DimVector<3> v12 = Vect<3>(V1, V2);
		DimVector<3> v13 = Vect<3>(V1, V3);
		DimVector<3> v23 = Vect<3>(V2, V3);
		_diameter = max({ v12.norm(), v13.norm(), v23.norm() });

		_measure = 0.5 * v12.cross(v13).norm();

		_center = DomPoint((V1->X + V2->X + V3->X) / 3, (V1->Y + V2->Y + V3->Y) / 3, (V1->Z + V2->Z + V3->Z) / 3);

		_detJacobian = _measure / RefShape()->Measure();

		/*DimMatrix<2> inverseJacobian;
		inverseJacobian(0, 0) = (V3->Y - V1->Y) / ((V3->Y - V1->Y)*(V2->X - V1->X) - (V3->X - V1->X)*(V2->Y - V1->Y));
		inverseJacobian(0, 1) = -(V3->X - V1->X) / ((V3->Y - V1->Y)*(V2->X - V1->X) - (V3->X - V1->X)*(V2->Y - V1->Y));
		inverseJacobian(1, 0) = (V2->Y - V1->Y) / ((V2->Y - V1->Y)*(V3->X - V1->X) - (V2->X - V1->X)*(V3->Y - V1->Y));
		inverseJacobian(1, 1) = -(V2->X - V1->X) / ((V2->Y - V1->Y)*(V3->X - V1->X) - (V2->X - V1->X)*(V3->Y - V1->Y));
		_inverseJacobianTranspose = inverseJacobian.transpose();*/
	}

	ReferenceShape<2>* RefShape() const
	{
		return &TriangleShape::RefTriangle;
	}
	
	inline vector<Vertex*> Vertices() const override
	{
		return vector<Vertex*> { V1, V2, V3 };
	}

	static ReferenceTriangle* InitReferenceShape()
	{
		return &TriangleShape::RefTriangle;
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
		return Geometry::IsInTriangle(*V1, *V2, *V3, p, _measure);
	}

	inline double DetJacobian() const
	{
		return _detJacobian;
	}
	inline DimMatrix<2> InverseJacobianTranspose() const
	{
		//return _inverseJacobianTranspose;
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
		/*double x = domainPoint.X;
		double y = domainPoint.Y;

		double t = ((V3->Y - V1->Y)*(x - V1->X) - (V3->X - V1->X)*(y - V1->Y)) / ((V3->Y - V1->Y)*(V2->X - V1->X) - (V3->X - V1->X)*(V2->Y - V1->Y));
		double u = ((V2->Y - V1->Y)*(x - V1->X) - (V2->X - V1->X)*(y - V1->Y)) / ((V2->Y - V1->Y)*(V3->X - V1->X) - (V2->X - V1->X)*(V3->Y - V1->Y));
		RefPoint p(t, u);
		return p;*/
		assert(false);
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

	}
};