#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/partition_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/mpq_class.h>
#include "../Mesh/Vertex.h"
using namespace std;

// CGAL types
//typedef CGAL::Exact_predicates_exact_constructions_kernel   exactKernel;
typedef CGAL::Simple_cartesian<mpq_class>                   exactKernel; // must use this kernel because the function intersection is not thread-safe (see https://github.com/CGAL/cgal/issues/2685)
typedef CGAL::Exact_predicates_inexact_constructions_kernel inexactKernel;

class CGALWrapper
{
public:
	template <typename CGALKernel>
	static CGAL::Polygon_2<CGALKernel> CreatePolygonWithColinearityChecking(const vector<Vertex*>& vertices)
	{
		RotatingList<Vertex*> vert(vertices);
		CGAL::Polygon_2<CGALKernel> cgalPolygon;
		for (int i = 0; i < vertices.size(); ++i)
		{
			Vertex* v = vert.Get();
			DimVector<2> v1 = Vect<2>(vert.GetPrevious(), v);
			DimVector<2> v2 = Vect<2>(v, vert.GetNext());
			if (!AreCollinear(v1, v2))
				cgalPolygon.push_back(CGAL::Point_2<CGALKernel>(v->X, v->Y));
			vert.MoveNext();
		}
		return cgalPolygon;
	}

	template <typename CGALKernel>
	static CGAL::Polygon_2<CGALKernel> CreatePolygon(const vector<Vertex*>& vertices, bool checkColinearity = false)
	{
		if (checkColinearity)
			return CreatePolygonWithColinearityChecking<CGALKernel>(vertices);

		CGAL::Polygon_2<CGALKernel> cgalPolygon;
		for (Vertex* v : vertices)
			cgalPolygon.push_back(CGAL::Point_2<CGALKernel>(v->X, v->Y));
		return cgalPolygon;
	}

	template <typename CGALKernel>
	static void FillPolygon(const vector<Vertex*>& vertices, CGAL::Polygon_2<CGALKernel>& polyToFill)
	{
		assert(vertices.size() >= 3);
		for (Vertex* v : vertices)
			polyToFill.push_back(CGAL::Point_2<CGALKernel>(v->X, v->Y));
	}

	// Get vertices from CGAL polygon
	template <typename CGALKernel>
	static vector<Vertex*> ToVertices(const CGAL::Polygon_2<CGALKernel>& poly)
	{
		vector<Vertex*> vertices;
		for (auto it = poly.vertices_begin(); it != poly.vertices_end(); it++)
		{
			CGAL::Point_2<CGALKernel> p = *it;
			Vertex* v = new Vertex(0, CGAL::to_double(p.x()), CGAL::to_double(p.y()));
			vertices.push_back(v);
		}
		return vertices;
	}

	static vector<Vertex*> ToVertices(const CGAL::Partition_traits_2<inexactKernel>::Polygon_2& poly)
	{
		vector<Vertex*> vertices;
		for (auto it = poly.vertices_begin(); it != poly.vertices_end(); it++)
		{
			CGAL::Point_2<inexactKernel> p = *it;
			Vertex* v = new Vertex(0, p.x(), p.y());
			vertices.push_back(v);
		}
		return vertices;
	}

	static vector<Vertex*> ToVertices(const CGAL::Partition_traits_2<exactKernel>::Polygon_2& poly)
	{
		vector<Vertex*> vertices;
		for (auto it = poly.vertices_begin(); it != poly.vertices_end(); it++)
		{
			CGAL::Point_2<exactKernel> p = *it;
			Vertex* v = new Vertex(0, CGAL::to_double(p.x()), CGAL::to_double(p.y()));
			vertices.push_back(v);
		}
		return vertices;
	}

	template <typename CGALKernel>
	static bool CGALPolyContains(const CGAL::Polygon_2<CGALKernel>& poly, const CGAL::Point_2<CGALKernel>& p)
	{
		switch (CGAL::bounded_side_2(poly.vertices_begin(), poly.vertices_end(), p))
		{
		case CGAL::ON_BOUNDED_SIDE: // inside
			return true;
		case CGAL::ON_BOUNDARY: // on the boundary
			return true;
		case CGAL::ON_UNBOUNDED_SIDE: // outside
			return false;
		}
		return false;
	}
};