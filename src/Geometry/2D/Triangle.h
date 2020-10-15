#pragma once
#include "../../Mesh/Vertex.h"
#include "ReferenceTriangle.h"
#include "../PhysicalShapeWithConstantJacobian.h"
using namespace std;

class Triangle : public PhysicalShapeWithConstantJacobian<2>
{
private:
	DomPoint v1;
	DomPoint v2;
	DomPoint v3;

	double _diameter;
	double _measure;
	DomPoint _center;
	double _inRadius;

	DimMatrix<2> _inverseJacobianTranspose;
	double _detJacobian;

public:
	static ReferenceTriangle RefTriangle;

	Triangle() {}

	Triangle(const DomPoint& p1, const DomPoint& p2, const DomPoint& p3)
		: v1(p1), v2(p2), v3(p3)
	{
		Init();
	}

	inline DomPoint V1() const { return v1; }
	inline DomPoint V2() const { return v2; }
	inline DomPoint V3() const { return v3; }

	Triangle(const Triangle& shape) = default;

	void Init()
	{
		double lengthEdge12 = sqrt(pow(v2.X - v1.X, 2) + pow(v2.Y - v1.Y, 2));
		double lengthEdge23 = sqrt(pow(v3.X - v2.X, 2) + pow(v3.Y - v2.Y, 2));
		double lengthEdge13 = sqrt(pow(v3.X - v1.X, 2) + pow(v3.Y - v1.Y, 2));
		_diameter = max(lengthEdge12, max(lengthEdge23, lengthEdge13));

		_measure = 0.5 * abs(v1.X * (v2.Y - v3.Y) + v2.X * (v3.Y - v1.Y) + v3.X * (v1.Y - v2.Y));

		_inRadius = 2 * _measure / (lengthEdge12 + lengthEdge23 + lengthEdge13);

		_center = DomPoint((v1.X + v2.X + v3.X) / 3, (v1.Y + v2.Y + v3.Y) / 3);

		_detJacobian = _measure / RefTriangle.Measure();

		DimMatrix<2> inverseJacobian;
		inverseJacobian(0, 0) = (v3.Y - v1.Y) / ((v3.Y - v1.Y)*(v2.X - v1.X) - (v3.X - v1.X)*(v2.Y - v1.Y));
		inverseJacobian(0, 1) = -(v3.X - v1.X) / ((v3.Y - v1.Y)*(v2.X - v1.X) - (v3.X - v1.X)*(v2.Y - v1.Y));
		inverseJacobian(1, 0) = (v2.Y - v1.Y) / ((v2.Y - v1.Y)*(v3.X - v1.X) - (v2.X - v1.X)*(v3.Y - v1.Y));
		inverseJacobian(1, 1) = -(v2.X - v1.X) / ((v2.Y - v1.Y)*(v3.X - v1.X) - (v2.X - v1.X)*(v3.Y - v1.Y));
		_inverseJacobianTranspose = inverseJacobian.transpose();
	}

	PhysicalShape<2>* CreateCopy() const
	{
		return new Triangle(*this);
	}

	ReferenceShape<2>* RefShape() const
	{
		return &RefTriangle;
	}
	
	inline vector<DomPoint> Vertices() const override
	{
		return vector<DomPoint>{ v1, v2, v3 };
	}

	bool IsDegenerated() const override
	{
		DimVector<2> V1V2 = Vect<2>(v1, v2);
		DimVector<2> V1V3 = Vect<2>(v1, v3);
		double cross = V1V2[0] * V1V3[1] - V1V2[1] * V1V3[0];
		return abs(cross) < Point::Tolerance;
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
		return TriangleContains(v1, v2, v3, p, _measure);
	}

	static bool TriangleContains(const DomPoint& A, const DomPoint& B, const DomPoint& C, const DomPoint& P, double triangleArea)
	{
		// From https://math.stackexchange.com/questions/4322/check-whether-a-point-is-within-a-3d-triangle

		DimVector<3> PA = Vect<3>(P, A);
		DimVector<3> PB = Vect<3>(P, B);
		DimVector<3> PC = Vect<3>(P, C);
		// Barycentric coordinates
		double alpha = PB.cross(PC).norm() / (2 * triangleArea);
		double beta = PC.cross(PA).norm() / (2 * triangleArea);
		double gamma = PA.cross(PB).norm() / (2 * triangleArea);

		double tol = Utils::Eps;
		return alpha + tol > 0 && alpha < 1 + tol    // alpha >= 0 && alpha <= 1
			&& beta + tol > 0 && beta < 1 + tol    // beta  >= 0 && beta  <= 1
			&& gamma + tol > 0 && gamma < 1 + tol    // gamma >= 0 && gamma <= 1
			&& (abs(alpha + beta + gamma - 1) < tol); // alpha + beta + gamma = 1
	}

	inline double DetJacobian() const
	{
		return _detJacobian;
	}
	inline DimMatrix<2> InverseJacobianTranspose() const
	{
		return _inverseJacobianTranspose;
	}

	DomPoint ConvertToDomain(const RefPoint& refPoint) const
	{
		double t = refPoint.X;
		double u = refPoint.Y;

		DomPoint p;
		p.X = (v2.X - v1.X) * t + (v3.X - v1.X)*u + v1.X;
		p.Y = (v2.Y - v1.Y) * t + (v3.Y - v1.Y)*u + v1.Y;
		return p;
	}

	RefPoint ConvertToReference(const DomPoint& domainPoint) const
	{
		/*Vertex* V1 = _vertices[0];
		Vertex* V2 = _vertices[1];
		Vertex* V3 = _vertices[2];*/

		double x = domainPoint.X;
		double y = domainPoint.Y;

		double t = ((v3.Y - v1.Y)*(x - v1.X) - (v3.X - v1.X)*(y - v1.Y)) / ((v3.Y - v1.Y)*(v2.X - v1.X) - (v3.X - v1.X)*(v2.Y - v1.Y));
		double u = ((v2.Y - v1.Y)*(x - v1.X) - (v2.X - v1.X)*(y - v1.Y)) / ((v2.Y - v1.Y)*(v3.X - v1.X) - (v2.X - v1.X)*(v3.Y - v1.Y));
		RefPoint p(t, u);
		return p;
	}

	void ExportToMatlab(string color = "r") const override
	{
		MatlabScript script;
		script.PlotTriangle(v1, v2, v3, color);
	}

	void Serialize(ostream& os) const override
	{
		os << "Triangle";
		os << " ";
		v1.Serialize(os, 2);
		os << "--";
		v2.Serialize(os, 2);
		os << "--";
		v3.Serialize(os, 2);
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
		DomPoint lowerLeft(-1, -1);
		DomPoint lowerRight(1, -1);
		DomPoint upperLeft(-1, 1);

		Triangle t(lowerLeft, lowerRight, upperLeft);

		t.UnitTests();

		RefPoint llRef = t.ConvertToReference(lowerLeft);
		assert(llRef == RefPoint(0, 0));
		DomPoint llDom = t.ConvertToDomain(RefPoint(0, 0));
		assert(lowerLeft == llDom);

		RefPoint ulRef = t.ConvertToReference(upperLeft);
		assert(ulRef == RefPoint(0, 1));
	}
};

ReferenceTriangle Triangle::RefTriangle = ReferenceTriangle();