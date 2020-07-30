#pragma once
#include "../Mesh/2D/QuadrilateralShape.h"
using namespace std;

class Geometry
{
public:

	static QuadrilateralShape* CreateBoundingBox(vector<Vertex*> vertices)
	{
		double maxX = -INFINITY;
		double maxY = -INFINITY;
		double minX = INFINITY;
		double minY = INFINITY;
		for (Vertex* v : vertices)
		{
			if (v->X > maxX)
				maxX = v->X;
			if (v->Y > maxY)
				maxY = v->Y;
			if (v->X < minX)
				minX = v->X;
			if (v->Y < minY)
				minY = v->Y;
		}

		int number = -1;
		Vertex* lowerLeft = new Vertex(number, minX, minY);
		Vertex* lowerRight = new Vertex(number, maxX, minY);
		Vertex* upperRight = new Vertex(number, maxX, maxY);
		Vertex* upperLeft = new Vertex(number, minX, maxY);

		QuadrilateralShape* boundingRectangle = new QuadrilateralShape(lowerLeft, lowerRight, upperRight, upperLeft);
		return boundingRectangle;
	}

	template <int Dim>
	bool AreInDirectOrder(Vertex* A, Vertex* B, Vertex* C)
	{
		if (Dim == 2)
		{
			DimVector<2> AB = *B - *A;
			DimVector<2> AC = *C - *A;
			double angle = atan2(AC[1], AC[0]) - atan2(AB[1], AB[0]);
			if (angle < 0) // to get angle in [0, 2*Pi]
				angle += 2 * M_PI;

			return angle < M_PI;
		}
		else
			assert(false);
	}

	static bool IsInSegment(DomPoint A, DomPoint B, DomPoint P)
	{
		DimVector<2> AB = B - A;
		DimVector<2> AP = P - A;
		double AB_dot_AP = AB.dot(AP);
		return AB_dot_AP > 0 && AB_dot_AP < AB.dot(AB) && abs(AB_dot_AP - AB.norm()*AP.norm()) < Point::Tolerance;
	}

	static bool IsInTriangle(DomPoint A, DomPoint B, DomPoint C, DomPoint P, double triangleArea, double triangleDiameter)
	{
		// From https://math.stackexchange.com/questions/4322/check-whether-a-point-is-within-a-3d-triangle

		DimVector<3> PA = Vect<3>(P, A);
		DimVector<3> PB = Vect<3>(P, B);
		DimVector<3> PC = Vect<3>(P, C);
		// Barycentric coordinates
		double alpha = PB.cross(PC).norm() / (2 * triangleArea);
		double beta = PC.cross(PA).norm() / (2 * triangleArea);
		double gamma = PA.cross(PB).norm() / (2 * triangleArea);

		double tol = Point::Tolerance / triangleDiameter;
		return alpha + tol > 0 && alpha < 1 + tol    // alpha >= 0 && alpha <= 1
			&& beta  + tol > 0 && beta  < 1 + tol    // beta  >= 0 && beta  <= 1
			&& gamma + tol > 0 && gamma < 1 + tol    // gamma >= 0 && gamma <= 1
			&& (abs(alpha + beta + gamma -1) < tol); // alpha + beta + gamma = 1
	}

	static bool IsInTetrahedron(DomPoint A, DomPoint B, DomPoint C, DomPoint D, DomPoint P, double tetraVolume, double tetraDiameter)
	{
		// Barycentric coordinates
		double alpha = VolumeTetra(P, B, C, D) / tetraVolume;
		double beta = VolumeTetra(A, P, C, D) / tetraVolume;
		double gamma = VolumeTetra(A, B, P, D) / tetraVolume;
		double delta = VolumeTetra(A, B, C, P) / tetraVolume;

		double tol = Point::Tolerance / tetraDiameter;
		return alpha + tol > 0 && alpha < 1 + tol    // alpha >= 0 && alpha <= 1
			&& beta  + tol > 0 && beta  < 1 + tol    // beta  >= 0 && beta  <= 1
			&& gamma + tol > 0 && gamma < 1 + tol    // gamma >= 0 && gamma <= 1
			&& delta + tol > 0 && delta < 1 + tol    // delta >= 0 && delta <= 1
			&& (abs(alpha + beta + gamma + delta - 1) < tol); // alpha + beta + gamma + delta = 1
	}

	static double VolumeTetra(DomPoint A, DomPoint B, DomPoint C, DomPoint D)
	{
		DimMatrix<3> m;
		m.col(0) = Vect<3>(A, B);
		m.col(1) = Vect<3>(A, C);
		m.col(2) = Vect<3>(A, D);
		return abs(m.determinant()) / 6;
	}

};