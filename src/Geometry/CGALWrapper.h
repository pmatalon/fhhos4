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
private:
	static void my_cgal_failure_handler(
		const char *type,
		const char *expr,
		const char* file,
		int line,
		const char* msg)
	{
		Utils::Error("CGAL error!");
	}

public:
	static void Configure()
	{
		/*CGAL::Failure_function prev; // save the CGAL failure function
		prev = CGAL::set_error_handler(my_cgal_failure_handler); // replace it with my own
		CGAL::set_error_handler(prev); // put the old one back*/
		CGAL::set_error_handler(my_cgal_failure_handler);
	}

	template <typename CGALKernel>
	static CGAL::Polygon_2<CGALKernel> CreatePolygonWithColinearityChecking(const vector<DomPoint>& vertices)
	{
		RotatingList<DomPoint> vert(vertices);
		CGAL::Polygon_2<CGALKernel> cgalPolygon;
		for (int i = 0; i < vertices.size(); ++i)
		{
			DomPoint v = vert.Get();
			DimVector<2> v1 = Vect<2>(vert.GetPrevious(), v);
			DimVector<2> v2 = Vect<2>(v, vert.GetNext());
			if (!AreCollinear(v1, v2))
				cgalPolygon.push_back(CGAL::Point_2<CGALKernel>(v.X, v.Y));
			vert.MoveNext();
		}
		return cgalPolygon;
	}

	template <typename CGALKernel>
	static CGAL::Polygon_2<CGALKernel> CreatePolygon(const vector<DomPoint>& vertices, bool checkColinearity = false)
	{
		if (checkColinearity)
			return CreatePolygonWithColinearityChecking<CGALKernel>(vertices);

		CGAL::Polygon_2<CGALKernel> cgalPolygon;
		for (const DomPoint& v : vertices)
			cgalPolygon.push_back(CGAL::Point_2<CGALKernel>(v.X, v.Y));
		return cgalPolygon;
	}

	template <typename CGALKernel>
	static void FillPolygon(const vector<DomPoint>& vertices, CGAL::Polygon_2<CGALKernel>& polyToFill)
	{
		assert(vertices.size() >= 3);
		for (const DomPoint& v : vertices)
			polyToFill.push_back(CGAL::Point_2<CGALKernel>(v.X, v.Y));
	}

	// Get vertices from CGAL polygon
	template <typename CGALKernel>
	static vector<DomPoint> ToVertices(const CGAL::Polygon_2<CGALKernel>& poly)
	{
		vector<DomPoint> vertices;
		for (auto it = poly.vertices_begin(); it != poly.vertices_end(); it++)
		{
			CGAL::Point_2<CGALKernel> p = *it;
			vertices.emplace_back(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
		}
		return vertices;
	}

	static vector<DomPoint> ToVertices(const CGAL::Partition_traits_2<inexactKernel>::Polygon_2& poly)
	{
		vector<DomPoint> vertices;
		for (auto it = poly.vertices_begin(); it != poly.vertices_end(); it++)
		{
			CGAL::Point_2<inexactKernel> p = *it;
			vertices.emplace_back(p.x(), p.y());
		}
		return vertices;
	}

	static vector<DomPoint> ToVertices(const CGAL::Partition_traits_2<exactKernel>::Polygon_2& poly)
	{
		vector<DomPoint> vertices;
		for (auto it = poly.vertices_begin(); it != poly.vertices_end(); it++)
		{
			CGAL::Point_2<exactKernel> p = *it;
			vertices.emplace_back(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
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

	template <typename CGALKernel>
	static vector<CGAL::Polygon_2<CGALKernel>> ToSimplePolygons(const CGAL::Polygon_2<CGALKernel>& poly)
	{
		vector<CGAL::Polygon_2<CGALKernel>> simplePolys;

		vector<CGAL::Point_2<CGALKernel>> points;
		CGAL::Point_2<CGALKernel> doublon;
		for (auto it = poly.vertices_begin(); it != poly.vertices_end(); it++)
		{
			CGAL::Point_2<CGALKernel> p = *it;
			if (find(points.begin(), points.end(), p) != points.end())
			{
				doublon = p;
				break;
			}
			points.push_back(p);
		}
		typename CGAL::Polygon_2<CGALKernel>::Vertex_const_circulator vertex = poly.vertices_circulator();
		while (*vertex != doublon)
			vertex++;

		CGAL::Polygon_2<CGALKernel> simplePoly1;
		do {
			simplePoly1.push_back(*vertex);
			vertex++;
		} while (*vertex != doublon);

		CGAL::Polygon_2<CGALKernel> simplePoly2;
		do {
			simplePoly2.push_back(*vertex);
			vertex++;
		} while (*vertex != doublon);

		if (simplePoly1.size() > 2)
		{
			if (simplePoly1.is_simple())
				simplePolys.push_back(simplePoly1);
			else
				simplePolys = ToSimplePolygons(simplePoly1);
		}

		if (simplePoly2.size() > 2)
		{
			if (simplePoly2.is_simple())
				simplePolys.push_back(simplePoly2);
			else
			{
				auto otherSimplePolys = ToSimplePolygons(simplePoly2);
				for (auto sp : otherSimplePolys)
					simplePolys.push_back(sp);
			}
		}

		return simplePolys;
	}
};