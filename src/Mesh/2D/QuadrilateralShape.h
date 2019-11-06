#pragma once
#include "../Vertex.h"
#include "../ReferenceCartesianShape.h"
#include "../GeometricShapeWithReferenceShape.h"
using namespace std;

class QuadrilateralShape : public GeometricShapeWithReferenceShape<2>
{
private:
	double _diameter;
	double _measure;
	DomPoint _center;

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

	QuadrilateralShape(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4)
	{
		V1 = v1;
		V2 = v2;
		V3 = v3;
		V4 = v4;
		Init();
	}

	void Init()
	{
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

		double diag13 = (*V3 - *V1).norm();
		double diag24 = (*V4 - *V2).norm();
		double edge12 = (*V2 - *V1).norm();
		double edge23 = (*V3 - *V2).norm();
		double edge34 = (*V4 - *V3).norm();
		double edge41 = (*V1 - *V4).norm();

		_diameter = max(diag13, diag24);

		_measure = 0.25 * sqrt(diag13*diag13 * diag24*diag24 - pow(edge12*edge12 + edge34*edge34 - edge23*edge23 - edge41*edge41, 2));

		_center = DomPoint((V1->X + V2->X + V3->X + V4->X) / 4, (V1->Y + V2->Y + V3->Y + V4->Y) / 4);
	}

	ReferenceShape<2>* RefShape() const
	{
		return &RefSquare;
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

	double DetJacobian() const
	{
		assert(false);
	}
	DimMatrix<2> InverseJacobianTranspose() const
	{
		assert(false);
	}
	int RefFaceNumber(const GeometricShape<1>* face) const
	{
		assert(false);
	}

	double DetJacobian(RefPoint p) const
	{
		DimMatrix<2> jacobianMatrix = JacobianMatrix(p);
		return jacobianMatrix.determinant();
	}
	DimMatrix<2> InverseJacobianTranspose(RefPoint p) const
	{
		DimMatrix<2> jacobianMatrix = JacobianMatrix(p);
		return jacobianMatrix.inverse().transpose();
	}
	DimMatrix<2> JacobianMatrix(RefPoint p) const
	{
		double t = p.X;
		double u = p.Y;
		DimMatrix<2> jacobianMatrix;
		jacobianMatrix(0, 0) = 0.25 * (V1->X *(u-1) + V2->X * (1-u)  + V3->X*(u+1) + V4->X*(-u-1)); // sum { x_i * d(Khi_i)/dt }
		jacobianMatrix(0, 1) = 0.25 * (V1->X *(t-1) + V2->X * (-t-1) + V3->X*(t+1) + V4->X*(1-t)); // sum { x_i * d(Khi_i)/du }
		jacobianMatrix(1, 0) = 0.25 * (V1->Y *(u-1) + V2->Y * (1-u)  + V3->Y*(u+1) + V4->Y*(-u-1)); // sum { y_i * d(Khi_i)/dt }
		jacobianMatrix(1, 1) = 0.25 * (V1->Y *(t-1) + V2->Y * (-t-1) + V3->Y*(t+1) + V4->Y*(1-t)); // sum { y_i * d(Khi_i)/du }
		return jacobianMatrix;
	}



	// Formulas in the book by Elman et al. "Finite Elements and Fast Iterative Solvers"
	DomPoint ConvertToDomain(RefPoint refPoint) const
	{
		double t = refPoint.X;
		double u = refPoint.Y;

		DomPoint p;
		p.X = V1->X * Khi1(t, u) + V2->X * Khi2(t, u) + V3->X * Khi3(t, u) + V4->X * Khi4(t, u);
		p.Y = V1->Y * Khi1(t, u) + V2->Y * Khi2(t, u) + V3->Y * Khi3(t, u) + V4->Y * Khi4(t, u);
		return p;
	}

	RefPoint ConvertToReference(DomPoint domainPoint) const
	{
		//assert(false);
		// x = a0 + a1*t + a2*u + a3*t*u
		// y = b0 + b1*t + b2*u + b3*t*u
		double t, u;

		double x0 = -(a0 - domainPoint.X);
		double y0 = -(b0 - domainPoint.Y);

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
				if (u < -1 || u > 1)
					u = (-B - sqrt(discr)) / (2 * A);
			}
			else if (discr == 0)
				u = -B / (2 * A);
			else
				assert(false);

			if (a1 + a3 * u != 0)
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
		return RefPoint(t, u);
	}

private:
	inline static double Khi1(double t, double u)
	{
		return (t - 1)*(u - 1) / 4;
	}
	inline static double Khi2(double t, double u)
	{
		return -(t + 1)*(u - 1) / 4;
	}
	inline static double Khi3(double t, double u)
	{
		return (t + 1)*(u + 1) / 4;
	}
	inline static double Khi4(double t, double u)
	{
		return -(t - 1)*(u + 1) / 4;
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

	static void UnitTests()
	{
		int number = 0;
		Vertex lowerLeft(number, -1, -1);
		Vertex lowerRight(number, 1, -1);
		Vertex upperRight(number, 2, 2);
		Vertex upperLeft(number, -1, 1);

		QuadrilateralShape q(&lowerLeft, &lowerRight, &upperRight, &upperLeft);
		RefPoint llRef = q.ConvertToReference(lowerLeft);
		assert(llRef == RefPoint(-1, -1));
		DomPoint llDom = q.ConvertToDomain(RefPoint(-1, -1));
		assert(lowerLeft == llDom);

		RefPoint urRef = q.ConvertToReference(upperRight);
		assert(urRef == RefPoint(1, 1));
	}
};

ReferenceCartesianShape<2> QuadrilateralShape::RefSquare = ReferenceCartesianShape<2>();