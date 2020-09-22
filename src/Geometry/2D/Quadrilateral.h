#pragma once
#include "../../Mesh/Vertex.h"
#include "../ReferenceCartesianShape.h"
#include "../PhysicalShape.h"
using namespace std;

class Quadrilateral : public PhysicalShape<2>
{
private:
	double _diameter;
	double _measure;
	DomPoint _center;
	double _inRadius;

	double a0;
	double a1;
	double a2;
	double a3;

	double b0;
	double b1;
	double b2;
	double b3;

public:
	Vertex* V1;
	Vertex* V2;
	Vertex* V3;
	Vertex* V4;

	static ReferenceCartesianShape<2> RefSquare;

	Quadrilateral(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4)
	{
		V1 = v1;
		V2 = v2;
		V3 = v3;
		V4 = v4;
		Init();
	}

	// Copy constructor
	Quadrilateral(const Quadrilateral& shape) = default;

	void Init()
	{
		double diag13 = (*V3 - *V1).norm();
		double diag24 = (*V4 - *V2).norm();
		double edge12 = (*V2 - *V1).norm();
		double edge23 = (*V3 - *V2).norm();
		double edge34 = (*V4 - *V3).norm();
		double edge41 = (*V1 - *V4).norm();

		_diameter = max(diag13, diag24);

		_measure = 0.25 * sqrt(4 * diag13*diag13 * diag24*diag24 - pow(edge12*edge12 + edge34 * edge34 - edge23 * edge23 - edge41 * edge41, 2));

		_inRadius = 2 * _measure / (edge12 + edge23 + edge34 + edge41);

		_center = DomPoint((V1->X + V2->X + V3->X + V4->X) / 4, (V1->Y + V2->Y + V3->Y + V4->Y) / 4);

		// x = a0 + a1*t + a2*u + a3*t*u
		// y = b0 + b1*t + b2*u + b3*t*u
		a0 = 0.25 * ( V1->X + V2->X + V3->X + V4->X);
		a1 = 0.25 * (-V1->X + V2->X + V3->X - V4->X);
		a2 = 0.25 * (-V1->X - V2->X + V3->X + V4->X);
		a3 = 0.25 * ( V1->X - V2->X + V3->X - V4->X);

		b0 = 0.25 * ( V1->Y + V2->Y + V3->Y + V4->Y);
		b1 = 0.25 * (-V1->Y + V2->Y + V3->Y - V4->Y);
		b2 = 0.25 * (-V1->Y - V2->Y + V3->Y + V4->Y);
		b3 = 0.25 * ( V1->Y - V2->Y + V3->Y - V4->Y);

		// Sometimes a3 or b3 can be really small, but not zero. Then we set it to zero be tested in ConvertToReference()
		if (abs(a3 / _measure) < 1e-10)
			a3 = 0;
		if (abs(b3 / _measure) < 1e-10)
			b3 = 0;
	}

	PhysicalShape<2>* CreateCopy() const
	{
		return new Quadrilateral(*this);
	}

	ReferenceShape<2>* RefShape() const
	{
		return &RefSquare;
	}

	inline vector<Vertex*> Vertices() const override
	{
		return vector<Vertex*> { V1, V2, V3, V4 };
	}

	bool IsDegenerated() const override
	{
		assert(false && "To implement");
	}

	void ReshapeByMovingIntersection(Vertex* oldIntersect, Vertex* newIntersect) override
	{
		if (*V1 == *oldIntersect)
			V1 = newIntersect;
		else if (*V2 == *oldIntersect)
			V2 = newIntersect;
		else if (*V3 == *oldIntersect)
			V3 = newIntersect;
		else if (*V4 == *oldIntersect)
			V4 = newIntersect;
		else
			assert(false && "This quadrilateral does not have this vertex.");
		Init();
	}

	static ReferenceCartesianShape<2>* InitReferenceShape()
	{
		return &RefSquare;
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
	inline bool Contains(const DomPoint& p) const override
	{
		assert(false && "Not implemented");
	}

	double DetJacobian(const RefPoint& p) const override
	{
		DimMatrix<2> jacobianMatrix = JacobianMatrix(p);
		return jacobianMatrix.determinant();
	}
	int DetJacobianDegree() const override
	{
		return 2;
	}
	DimMatrix<2> InverseJacobianTranspose(const RefPoint& p) const
	{
		DimMatrix<2> jacobianMatrix = JacobianMatrix(p);
		return jacobianMatrix.inverse().transpose();
	}
	DimMatrix<2> JacobianMatrix(const RefPoint& p) const
	{
		double t = p.X;
		double u = p.Y;

		DimMatrix<2> jacobianMatrix;
		jacobianMatrix(0, 0) = a1 + a3 * u;
		jacobianMatrix(0, 1) = a2 + a3 * t;
		jacobianMatrix(1, 0) = b1 + b3 * u;
		jacobianMatrix(1, 1) = b2 + b3 * t;

		return jacobianMatrix;
	}

	// Formulas in Silva et al. "Exact and efficient interpolation using finite elements shape functions" (2009)
	// where ksi = t, eta = u
	DomPoint ConvertToDomain(const RefPoint& refPoint) const
	{
		double t = refPoint.X;
		double u = refPoint.Y;

		DomPoint p;
		p.X = a0 + a1*t + a2*u + a3*t*u;
		p.Y = b0 + b1*t + b2*u + b3*t*u;
		return p;
	}

	RefPoint ConvertToReference(const DomPoint& domainPoint) const
	{
		double t, u;

		double x0 = domainPoint.X - a0;
		double y0 = domainPoint.Y - b0;

		if (a3 == 0 && b3 == 0)
		{
			t = (b2 * x0 - a2 * y0) / (b2 * a1 - a2 * b1);
			u = (y0 - b1 * t) / b2;
		}
		else
		{
			double A = a3 * b2 - a2 * b3;
			double B = x0 * b3 + a1 * b2 - (y0*a3 + a2 * b1);
			double C = x0 * b1 - y0 * a1;

			double discr = B * B - 4 * A*C;
			if (discr > 0)
			{
				u = (-B + sqrt(discr)) / (2 * A);
				//if (u < -1 || u > 1)
				if (abs(u) > 1.0001)
					u = (-B - sqrt(discr)) / (2 * A);
			}
			else if (discr == 0)
				u = -B / (2 * A);
			else
				assert(false);

			if (abs(a1 + a3 * u) > 1e-16) // if (a1 + a3 * u != 0)
				t = (x0 - a2 * u) / (a1 + a3 * u);
			else if (a1 != 0 && a3 != 0)
			{
				u = -a1 / a3;
				t = (y0*a3 + a1 * b2) / (a3 * b1 - a1 * b3);
			}
			else if (a1 == 0 && a3 == 0)
			{
				u = x0 / a2;
				t = (y0*a2 - x0 * b2) / (a2 * b1 + x0 * b3);
			}
			else
				assert(false);
		}
		//if (abs(t) > 1.1 || abs(u) > 1.1)
			//assert(false && "The point is not included in the quadrilateral");
		return RefPoint(t, u);
	}

public:

	void Serialize(ostream& os) const override
	{
		os << "Quadrilateral";
		os << " ";
		V1->Serialize(os, 2);
		os << "-";
		V2->Serialize(os, 2);
		os << "-";
		V3->Serialize(os, 2);
		os << "-";
		V4->Serialize(os, 2);
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	static void Test()
	{
		int number = 0;
		Vertex lowerLeft(number, -1, -1);
		Vertex lowerRight(number, 1, -1);
		Vertex upperRight(number, 2, 2);
		Vertex upperLeft(number, -1, 1);

		Quadrilateral q(&lowerLeft, &lowerRight, &upperRight, &upperLeft);

		q.UnitTests();

		RefPoint llRef = q.ConvertToReference(lowerLeft);
		assert(llRef == RefPoint(-1, -1));
		DomPoint llDom = q.ConvertToDomain(RefPoint(-1, -1));
		assert(lowerLeft == llDom);

		RefPoint urRef = q.ConvertToReference(upperRight);
		assert(urRef == RefPoint(1, 1));


		Vertex V1(number, 0, 0);
		Vertex V2(number, 0.13, 0);
		Vertex V3(number, 0.13, 0.05);
		Vertex V4(number, 0, 0.05);
		Quadrilateral q2(&V1, &V2, &V3, &V4);
		DomPoint dom = q.ConvertToDomain(RefPoint(-0.069222, 0.534611));
	}
};

ReferenceCartesianShape<2> Quadrilateral::RefSquare = ReferenceCartesianShape<2>();