#pragma once
#include "../Mesh/2D/TriangleShape.h"
#include "../Mesh/2D/QuadrilateralShape.h"
using namespace std;

class Geometry
{
public:
	static vector<TriangleShape*> Triangulation(vector<Vertex*> vertices)
	{
		vector<TriangleShape*> triangles;

		double sumX = 0;
		double sumY = 0;
		for (Vertex* v : vertices)
		{
			sumX += v->X;
			sumY += v->Y;
		}

		Vertex* center = new Vertex(0, sumX / vertices.size(), sumY / vertices.size());

		for (int i = 0; i < vertices.size(); i++)
		{
			TriangleShape* subTriangle = new TriangleShape(vertices[i], vertices[(i + 1) % vertices.size()], center);
			triangles.push_back(subTriangle);
		}
		return triangles;
	}

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

};