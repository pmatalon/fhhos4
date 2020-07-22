#pragma once
#include "../../Utils/Geometry.h"
#include "../CartesianShape.h"
#include "TriangleShape.h"
#include "../../Utils/RotatingList.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/partition_2.h>
#include <CGAL/Partition_traits_2.h>
#include <cassert>
#include <list>
using namespace std;

// CGAL types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Partition_traits_2<K>                         Traits;
typedef Traits::Polygon_2                                   Polygon_2;
typedef Traits::Point_2                                     Point_2;
typedef std::list<Polygon_2>                                Polygon_list;



class PolygonalShape : public GeometricShapeWithReferenceShape<2>
{
private:
	vector<Vertex*> _vertices;
	Polygon_2 _cgalPolygon; // Polygon of the CGAL library

	double _diameter;
	double _measure;
	Vertex* _center;
	double _inRadius;

	vector<GeometricShapeWithReferenceShape<2>*> _triangulation;
	QuadrilateralShape* _boundingBox;

public:

	PolygonalShape(vector<Vertex*> vertices) 
		: _vertices(vertices)
	{
		Init(vertices);
	}

	PolygonalShape(const PolygonalShape& shape) = default;

	void Init(vector<Vertex*> vertices)
	{
		assert(vertices.size() >= 3);
		_vertices = vertices;

		_cgalPolygon.clear();

		// Bug of CGAL: https://github.com/CGAL/cgal/issues/2575
		// The partitioning function used for the decomposition in convex sub-parts will fail if two edges are collinear.
		// So we need to remove the "useless" vertices for the creation of the CGAL polygon.
		RotatingList<Vertex*> vert(_vertices);
		for (int i = 0; i < _vertices.size(); ++i)
		{
			Vertex* v = vert.Get();
			DimVector<2> v1 = Vect<2>(vert.GetPrevious(), v);
			DimVector<2> v2 = Vect<2>(v, vert.GetNext());
			if (!AreCollinear(v1, v2))
				_cgalPolygon.push_back(Point_2(v->X, v->Y));
			vert.MoveNext();
		}
		//for (Vertex* v : vertices)
			//_cgalPolygon.push_back(Point_2(v->X, v->Y));
		assert(vertices.size() >= 3);

		Init();
	}

	void Init()
	{
		assert(_cgalPolygon.is_simple());
		assert(_cgalPolygon.is_counterclockwise_oriented());

		_diameter = 0;
		double sumX = 0;
		double sumY = 0;
		for (Vertex* v1 : _vertices)
		{
			for (Vertex* v2 : _vertices)
			{
				if (v1 != v2)
				{
					double diagonal12 = (*v2 - *v1).norm();
					if (diagonal12 > _diameter)
						_diameter = diagonal12;
				}
			}
			sumX += v1->X;
			sumY += v1->Y;
		}
		
		_center = new Vertex(0, sumX / _vertices.size(), sumY / _vertices.size());

		_triangulation.clear();
		Triangulation();

		_boundingBox = Geometry::CreateBoundingBox(_vertices);

		_measure = 0;
		for (GeometricShapeWithReferenceShape<2>* t : _triangulation)
			_measure += t->Measure();

		assert(abs(_measure - _cgalPolygon.area()) < 1e-12);

		// TODO
		_inRadius = 0;
	}

private:
	void Triangulation()
	{
		if (_vertices.size() == 3)
		{
			TriangleShape* triangle = new TriangleShape(_vertices[0], _vertices[1], _vertices[2]);
			_triangulation.push_back(triangle);
		}
		else if (_vertices.size() == 4 && this->IsConvex())
		{
			TriangleShape* triangle1 = new TriangleShape(_vertices[0], _vertices[1], _vertices[2]);
			_triangulation.push_back(triangle1);
			TriangleShape* triangle2 = new TriangleShape(_vertices[2], _vertices[3], _vertices[0]);
			_triangulation.push_back(triangle2);
		}
		else if (this->IsConvex())
			_triangulation = TriangulationByCenter(_vertices);
		else
			_triangulation = CGALTriangulation();
	}

	static vector<GeometricShapeWithReferenceShape<2>*> TriangulationByCenter(vector<Vertex*> vertices)
	{
		// Requirement: the polygon defined by the vertices must be convex!

		vector<GeometricShapeWithReferenceShape<2>*> triangles;

		double sumX = 0;
		double sumY = 0;
		for (Vertex* v : vertices)
		{
			sumX += v->X;
			sumY += v->Y;
		}

		Vertex* center = new Vertex(-1, sumX / vertices.size(), sumY / vertices.size());

		if (vertices.size() == 3)
		{
			TriangleShape* triangle = new TriangleShape(vertices[0], vertices[1], vertices[2]);
			triangles.push_back(triangle);
		}
		else
		{
			for (int i = 0; i < vertices.size(); i++)
			{
				TriangleShape* subTriangle = new TriangleShape(vertices[i], vertices[(i + 1) % vertices.size()], center);
				triangles.push_back(subTriangle);
			}
		}
		return triangles;
	}

	/*static void my_cgal_failure_handler(
		const char *type,
		const char *expr,
		const char* file,
		int line,
		const char* msg)
	{
		// report the error in some way.
	}*/

	vector<GeometricShapeWithReferenceShape<2>*> CGALTriangulation()
	{
		vector<GeometricShapeWithReferenceShape<2>*> triangulation;

		Polygon_list partition_polys;

		//CGAL::Failure_function prev; // save the CGAL failure function
		//prev = CGAL::set_error_handler(my_cgal_failure_handler); // replace it with my own
		CGAL::greene_approx_convex_partition_2(_cgalPolygon.vertices_begin(), _cgalPolygon.vertices_end(), back_inserter(partition_polys));
		//CGAL::set_error_handler(prev); // put the old one back

		for (Polygon_2 p : partition_polys)
		{
			vector<Vertex*> vertices = Vertices(p);
			if (p.size() == 3)
			{
				TriangleShape* subTriangle = new TriangleShape(vertices[0], vertices[1], vertices[2]);
				triangulation.push_back(subTriangle);
			}
			else
			{
				vector<GeometricShapeWithReferenceShape<2>*> subTriangles = TriangulationByCenter(vertices);
				for (auto tri : subTriangles)
					triangulation.push_back(tri);
			}
		}

		return triangulation;
	}

private:
	static vector<Vertex*> Vertices(Polygon_2& poly)
	{
		vector<Vertex*> vertices;
		for (auto it = poly.vertices_begin(); it != poly.vertices_end(); it++)
		{
			Point_2 p = *it;
			Vertex* v = new Vertex(0, p.x(), p.y());
			vertices.push_back(v);
		}
		return vertices;
	}

public:
	GeometricShapeWithReferenceShape<2>* CreateCopy() const
	{
		PolygonalShape* copy = new PolygonalShape(*this);
		//copy->_triangulation = Triangulation(copy->_vertices);
		//copy->_boundingBox = static_cast<QuadrilateralShape*>(this->_boundingBox->CreateCopy());
		copy->Init();
		return copy;
	}

	bool IsMadeOfSubShapes() const override
	{
		return true;
	}
	vector<GeometricShapeWithReferenceShape<2>*> SubShapes() const override
	{
		return _triangulation;
	}

	ReferenceShape<2>* RefShape() const
	{
		return &RectangleShape::RefCartShape;
	}

	inline vector<Vertex*> Vertices() const override
	{
		return _vertices;
	}

	bool IsDegenerated() const override
	{
		assert(false && "To implement");
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
		return *_center;
	}
	inline bool IsConvex() const override
	{
		/*if (_vertices.size() < 4)
			return true;
		// Not convex if a cord is outside the shape
		for (int i = 0; i < _vertices.size(); i++)
		{
			for (int j = i + 1; j < _vertices.size(); j++)
			{
				if (!Contains(DomPoint::Middle(*_vertices[i], *_vertices[j])))
					return false;
			}
		}
		return true;*/
		return _cgalPolygon.is_convex();
	}
	inline double InRadius() const override
	{
		return _inRadius;
	}
	inline bool Contains(DomPoint p) const override
	{
		/*for (GeometricShapeWithReferenceShape<2>* t : _triangulation)
		{
			if (t->Contains(p))
				return true;
		}
		return false;*/
		switch (CGAL::bounded_side_2(_cgalPolygon.vertices_begin(), _cgalPolygon.vertices_end(), Point_2(p.X, p.Y), K()))
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

	void ExportSubShapesToMatlab() const override
	{
		MatlabScript script;
		vector<string> options = { "r", "b", "g", "c", "m", "y" };
		for (int i = 0; i < _triangulation.size(); i++)
		{
			string option = options[i % options.size()];
			auto vertices = _triangulation[i]->Vertices();
			if (vertices.size() == 3)
				script.PlotTriangle(*vertices[0], *vertices[1], *vertices[2], option);
			else
				assert(false && "To implement");
		}
	}

	double Integral(RefFunction boundingBoxDefinedFunction) const override
	{
		// For the inner triangles, the bounding box is the domain
		DomFunction boundingBoxFunction = [this, boundingBoxDefinedFunction](DomPoint boundingBoxPoint)
		{
			RefPoint p = _boundingBox->ConvertToReference(boundingBoxPoint);
			return boundingBoxDefinedFunction(p);
		};

		double integral = 0;
		for (GeometricShapeWithReferenceShape<2>* t : _triangulation)
			integral += t->Integral(boundingBoxFunction);
		return integral;
	}

	double Integral(RefFunction boundingBoxDefinedFunction, int polynomialDegree) const override
	{
		// For the inner triangles, the bounding box is the domain
		DomFunction boundingBoxFunction = [this, boundingBoxDefinedFunction](DomPoint boundingBoxPoint)
		{
			RefPoint p = _boundingBox->ConvertToReference(boundingBoxPoint);
			return boundingBoxDefinedFunction(p);
		};

		double integral = 0;
		for (GeometricShapeWithReferenceShape<2>* t : _triangulation)
			integral += t->Integral(boundingBoxFunction, polynomialDegree);
		return integral;
	}

	inline double DetJacobian(RefPoint pointInReferenceSquare) const
	{
		return _boundingBox->DetJacobian(pointInReferenceSquare);
	}
	inline DimMatrix<2> InverseJacobianTranspose(RefPoint pointInReferenceSquare) const
	{
		return _boundingBox->InverseJacobianTranspose(pointInReferenceSquare);
	}
	inline int DetJacobianDegree() const
	{
		return _boundingBox->DetJacobianDegree();
	}

	DomPoint ConvertToDomain(RefPoint pointInReferenceSquare) const
	{
		return _boundingBox->ConvertToDomain(pointInReferenceSquare);
	}

	RefPoint ConvertToReference(DomPoint domainPoint) const
	{
		return _boundingBox->ConvertToReference(domainPoint);
	}

	void Serialize(ostream& os) const override
	{
		os << "Polygon";
	}

	~PolygonalShape()
	{
		for (GeometricShapeWithReferenceShape<2>* t : _triangulation)
			delete t;
		delete _center;
		delete _boundingBox;
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	static void Test()
	{
		int number = 0;
		Vertex lowerLeft(number, -1, -1);
		Vertex lowerRight(number, 1, -1);
		Vertex upperRight(number, 1, 1);
		Vertex upperLeft(number, -1, 1);
		vector<Vertex*> vertices{ &lowerLeft, &lowerRight, &upperRight, &upperLeft };
		PolygonalShape polygRefSquare(vertices);

		//--------------------------------------//
		polygRefSquare.UnitTests();

		//--------------------------------------//
		RefPoint llRef = polygRefSquare.ConvertToReference(lowerLeft);
		assert(llRef == RefPoint(-1, -1));
		DomPoint llDom = polygRefSquare.ConvertToDomain(RefPoint(-1, -1));
		assert(lowerLeft == llDom);

		//--------------------------------------//
		RefPoint urRef = polygRefSquare.ConvertToReference(upperRight);
		assert(urRef == RefPoint(1, 1));
		DomPoint urDom = polygRefSquare.ConvertToDomain(RefPoint(1, 1));
		assert(upperRight == urDom);

		//--------------------------------------//
		RefFunction refX = [](RefPoint p) { return p.X; };
		double integral = polygRefSquare.Integral(refX);
		assert(integral == 0);

		//--------------------------------------//
		RefFunction anyFunction = [](RefPoint p) { return p.X + p.Y*p.Y + 1; };

		RectangleShape realRefSquare(&lowerLeft, 2);
		double integralOverRealSquare = realRefSquare.Integral(anyFunction);
		double integralOverPolygonalSquare = polygRefSquare.Integral(anyFunction);
		assert(abs(integralOverRealSquare - integralOverPolygonalSquare) < 1e-10);

		//--------------------------------------//
		double h = 2;
		lowerLeft = Vertex(number, 0, 0);
		lowerRight = Vertex(number, h, 0);
		upperRight = Vertex(number, h, h);
		upperLeft = Vertex(number, 0, h);
		PolygonalShape polygSquare(vector<Vertex*>{ &lowerLeft, &lowerRight, &upperRight, &upperLeft });
		RectangleShape realSquare(&lowerLeft, h);
		integralOverRealSquare = realSquare.Integral(anyFunction);
		integralOverPolygonalSquare = polygSquare.Integral(anyFunction);
		assert(abs(integralOverRealSquare - integralOverPolygonalSquare) < 1e-10);

		//--------------------------------------//
		h = 1;
		lowerLeft = Vertex(number, 0, 0);
		lowerRight = Vertex(number, h, 0);
		upperRight = Vertex(number, h, h);
		upperLeft = Vertex(number, 0, h);
		PolygonalShape polygSquare1(vector<Vertex*>{ &lowerLeft, &lowerRight, &upperRight, &upperLeft });
		realSquare = RectangleShape(&lowerLeft, h);
		integralOverRealSquare = realSquare.Integral(anyFunction);
		integralOverPolygonalSquare = polygSquare1.Integral(anyFunction);
		assert(abs(integralOverRealSquare - integralOverPolygonalSquare) < 1e-10);
	}
};