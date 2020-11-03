#pragma once
#include "../../Utils/RotatingList.h"
#include "../../Utils/MatlabScript.h"
#include "../CartesianShape.h"
#include "../CGALWrapper.h"
#include "Triangle.h"
#include "Quadrilateral.h"
#include <CGAL/convex_hull_2.h>
#include <cassert>
#include <list>
#include <exception>
using namespace std;

typedef exactKernel chosenKernel;

class InvalidPolygonException : public exception
{
public:
	vector<DomPoint> Vertices;
	InvalidPolygonException(const vector<DomPoint>& vertices) : Vertices(vertices)
	{}
	virtual const char* what() const throw()
	{
		return "InvalidPolygonException";
	}
};

class Polygon : public PhysicalShape<2>
{
private:
	vector<DomPoint> _vertices;
	//vector<Vertex*> _nonColinearVertices;

	CGAL::Polygon_2<chosenKernel> _cgalPolygon; // Polygon of the CGAL library

	double _diameter;
	double _measure;
	DomPoint _center;
	double _inRadius;

	vector<Triangle> _triangulation;
	Quadrilateral _boundingBox;
	vector<DomPoint> _quadraturePoints;

public:
	Polygon() {}

	Polygon(const vector<DomPoint>& vertices, bool createTriangulationAndBoundingBox, bool checkForColinearVertices)
		: _vertices(vertices)
	{
		assert(vertices.size() >= 3);
		InitCGALPolygon(checkForColinearVertices);
		Init(createTriangulationAndBoundingBox);
	}

	Polygon(const Polygon& shape) = default;

	Polygon(const CGAL::Polygon_2<chosenKernel>& cgalPolygon, bool createTriangulationAndBoundingBox)
		: _cgalPolygon(cgalPolygon)
	{
		_vertices = CGALWrapper::ToVertices(cgalPolygon);
		Init(createTriangulationAndBoundingBox);
	}

	const CGAL::Polygon_2<chosenKernel>& CGALPolygon() const
	{
		return _cgalPolygon;
	}

private:
	void InitCGALPolygon(bool checkForColinearVertices)
	{
		// Bug of CGAL: https://github.com/CGAL/cgal/issues/2575
		// The partitioning function used for the decomposition in convex sub-parts will fail if two edges are collinear.
		// So we need to remove the "useless" vertices for the creation of the CGAL polygon.

		//_nonColinearVertices = checkForColinearVertices ? NonColinearVertices(_vertices) : _vertices;
		//CGALWrapper::FillPolygon(_nonColinearVertices, _cgalPolygon);
		CGALWrapper::FillPolygon(_vertices, _cgalPolygon);

		if (!_cgalPolygon.is_simple())
		{
			Utils::Error("Non simple polygon detected. Printed below for Matlab:");
			ExportToMatlab("m");
			throw new InvalidPolygonException(_vertices);
		}
		else if (!_cgalPolygon.is_counterclockwise_oriented())
		{
			Utils::Error("Clockwise-oriented polygon detected. Printed below for Matlab:");
			ExportToMatlab("m");
			throw new InvalidPolygonException(_vertices);
		}
	}

	static vector<Vertex*> NonColinearVertices(const vector<Vertex*>& vertices)
	{
		vector<Vertex*> nonColinearVertices;
		RotatingList<Vertex*> vert(vertices);
		for (int i = 0; i < vertices.size(); ++i)
		{
			Vertex* v = vert.Get();
			DimVector<2> v1 = Vect<2>(vert.GetPrevious(), v);
			DimVector<2> v2 = Vect<2>(v, vert.GetNext());
			if (!AreCollinear(v1, v2))
				nonColinearVertices.push_back(v);
			vert.MoveNext();
		}
		return nonColinearVertices;
	}

	void Init(bool createTriangulationAndBoundingBox)
	{
		_diameter = 0;
		double sumX = 0;
		double sumY = 0;
		for (const DomPoint& v1 : _vertices)
		{
			for (const DomPoint& v2 : _vertices)
			{
				if (v1 != v2)
				{
					double diagonal12 = (v2 - v1).norm();
					if (diagonal12 > _diameter)
						_diameter = diagonal12;
				}
			}
			sumX += v1.X;
			sumY += v1.Y;
		}
		
		_center = DomPoint(sumX / _vertices.size(), sumY / _vertices.size());

		_measure = CGAL::to_double(_cgalPolygon.area());

		if (createTriangulationAndBoundingBox)
		{
			ComputeTriangulation();
			ComputeBoundingBox();
		}

		// TODO
		_inRadius = 0;
	}

public:
	void ComputeTriangulation()
	{
		if (!_triangulation.empty())
			return;

		if (_vertices.size() <= 3 || this->IsConvex())
			_triangulation = ConvexTriangulation(_vertices);
		else
			_triangulation = NonConvexTriangulation();

		auto it = _triangulation.begin();
		while (it != _triangulation.end())
		{
			if (it->Measure() < Utils::Eps*this->_measure)
				it = _triangulation.erase(it);
			else
				it++;
		}

		assert(_triangulation.size() > 0);

		//InitQuadraturePoints();
	}

	void ComputeBoundingBox()
	{
		if (_boundingBox.Measure() > 0)
			return;
		_boundingBox = CreateBoundingBox(_vertices);
	}

	static Quadrilateral CreateBoundingBox(const vector<DomPoint>& vertices)
	{
		double maxX = -INFINITY;
		double maxY = -INFINITY;
		double minX = INFINITY;
		double minY = INFINITY;
		for (const DomPoint& v : vertices)
		{
			if (v.X > maxX)
				maxX = v.X;
			if (v.Y > maxY)
				maxY = v.Y;
			if (v.X < minX)
				minX = v.X;
			if (v.Y < minY)
				minY = v.Y;
		}

		DomPoint lowerLeft (minX, minY);
		DomPoint lowerRight(maxX, minY);
		DomPoint upperRight(maxX, maxY);
		DomPoint upperLeft (minX, maxY);

		Quadrilateral boundingRectangle(lowerLeft, lowerRight, upperRight, upperLeft);
		return boundingRectangle;
	}

	void Refine() override
	{
		if (!_triangulation.empty())
			return;
		ComputeTriangulation();
	}

private:
	static vector<Triangle> ConvexTriangulation(const vector<DomPoint>& vertices)
	{
		assert(vertices.size() > 2);
		vector<Triangle> triangles;
		if (vertices.size() == 3)
			triangles.emplace_back(vertices[0], vertices[1], vertices[2]);
		else if (vertices.size() == 4)
		{
			triangles.emplace_back(vertices[0], vertices[1], vertices[2]);
			triangles.emplace_back(vertices[2], vertices[3], vertices[0]);
		}
		else
		{
			triangles.emplace_back(vertices[0], vertices[1], vertices[2]);

			// remove vertices[1]
			vector<DomPoint> remainingVertices;
			for (int i = 2; i < vertices.size(); i++)
				remainingVertices.push_back(vertices[i]);
			remainingVertices.push_back(vertices[0]);

			vector<Triangle> otherTriangles = ConvexTriangulation(remainingVertices);
			for (auto t : otherTriangles)
				triangles.push_back(t);
		}
		return triangles;
	}
	static vector<PhysicalShape<2>*> BarycentricTriangulation(const vector<DomPoint>& vertices)
	{
		// Requirement: the polygon defined by the vertices must be convex!

		vector<PhysicalShape<2>*> triangles;

		double sumX = 0;
		double sumY = 0;
		for (const DomPoint& v : vertices)
		{
			sumX += v.X;
			sumY += v.Y;
		}

		DomPoint center(sumX / vertices.size(), sumY / vertices.size());

		if (vertices.size() == 3)
		{
			Triangle* triangle = new Triangle(vertices[0], vertices[1], vertices[2]);
			triangles.push_back(triangle);
		}
		else
		{
			for (int i = 0; i < vertices.size(); i++)
			{
				Triangle* subTriangle = new Triangle(vertices[i], vertices[(i + 1) % vertices.size()], center);
				triangles.push_back(subTriangle);
			}
		}
		return triangles;
	}

	vector<Triangle> NonConvexTriangulation()
	{
		vector<Triangle> triangulation;

		list<CGAL::Partition_traits_2<chosenKernel>::Polygon_2> convexPartition;

		try
		{
			//CGAL::greene_approx_convex_partition_2(_cgalPolygon.vertices_begin(), _cgalPolygon.vertices_end(), back_inserter(convexPartition));
			CGAL::optimal_convex_partition_2(_cgalPolygon.vertices_begin(), _cgalPolygon.vertices_end(), back_inserter(convexPartition));
		}
		catch (CGAL::Failure_exception e)
		{
			cout << endl;
			this->ExportCGALPolyToMatlab();
			Utils::FatalError("CGAL failed to compute a partitioning of the above polygon: " + e.message());
		}

		assert(convexPartition.size() > 1);

		for (CGAL::Partition_traits_2<chosenKernel>::Polygon_2 p : convexPartition)
		{
			vector<DomPoint> vertices = CGALWrapper::ToVertices(p);
			vector<Triangle> subTriangles = ConvexTriangulation(vertices);
			for (auto tri : subTriangles)
				triangulation.push_back(tri);
		}

		return triangulation;
	}

	void InitQuadraturePoints()
	{
		for (Triangle& t : _triangulation)
		{
			for (const RefPoint& refPoint : t.RefShape()->QuadraturePoints())
			{
				DomPoint domPoint = t.ConvertToDomain(refPoint);//this->ConvertToDomainAndSaveResult(refPoint);
				_quadraturePoints.push_back(domPoint);
			}
		}
	}

public:
	bool ConvexHullEmbeds(const PhysicalShape<2>* s) const override
	{
		vector<CGAL::Point_2<chosenKernel>> convexHullPoints;
		CGAL::convex_hull_2(_cgalPolygon.vertices_begin(), _cgalPolygon.vertices_end(), back_inserter(convexHullPoints));

		CGAL::Polygon_2<chosenKernel> convexHull(convexHullPoints.begin(), convexHullPoints.end());

		for (const DomPoint& v : s->Vertices())
		{
			if (!CGALWrapper::CGALPolyContains(convexHull, CGAL::Point_2<chosenKernel>(v.X, v.Y)))
				return false;
		}
		return true;
	}

	/*bool IntersectsWith(PhysicalShape<2>* other) override
	{
		Polygon* otherPolygon = dynamic_cast<Polygon*>(other);
		if (otherPolygon)
		{
			bool res = CGAL::do_intersect(_cgalPolygon, otherPolygon->_cgalPolygon);
			return res;
		}
		else
		{
			auto otherPolygon = CGALWrapper::CreatePolygon<exactKernel>(other->Vertices());
			bool res = CGAL::do_intersect(_cgalPolygon, otherPolygon);
			return res;
		}
	}*/

	
public:
	virtual PhysicalShape<2>* CreateCopy() const override
	{
		Polygon* copy = new Polygon(*this);
		return copy;
	}

	bool IsMadeOfSubShapes() const override
	{
		return true;
	}
	vector<const PhysicalShape<2>*> SubShapes() const override
	{
		assert(_triangulation.size() > 0);
		vector<const PhysicalShape<2>*> subShapes;
		for (const Triangle& t : _triangulation)
		{
			const PhysicalShape<2>* ps = &t;
			subShapes.push_back(ps);
		}
		return subShapes;
	}
	vector<PhysicalShape<2>*> SubShapes() override
	{
		assert(_triangulation.size() > 0);
		vector<PhysicalShape<2>*> subShapes;
		for (Triangle& t : _triangulation)
		{
			PhysicalShape<2>* ps = &t;
			subShapes.push_back(ps);
		}
		return subShapes;
	}

	const vector<Triangle>& Triangulation() const
	{
		return _triangulation;
	}

	virtual ReferenceShape<2>* RefShape() const override
	{
		return &RectangleShape::RefCartShape;
	}

	inline vector<DomPoint> Vertices() const override
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
		return _center;
	}
	inline DomPoint InteriorPoint() const override
	{
		if (IsConvex() || Contains(_center))
			return _center;
		assert(!_triangulation.empty());
		return _triangulation[0].Center();
	}
	inline bool IsConvex() const override
	{
		return _cgalPolygon.is_convex();
	}
	inline double InRadius() const override
	{
		return _inRadius;
	}
	inline bool Contains(const DomPoint& p) const override
	{
		return CGALWrapper::CGALPolyContains(_cgalPolygon, CGAL::Point_2<chosenKernel>(p.X, p.Y));
	}

public:
	void ExportSubShapesToMatlab() const override
	{
		assert(_triangulation.size() > 0);
		MatlabScript script;
		vector<string> options = { "r", "b", "g", "c", "m", "y" };
		for (int i = 0; i < _triangulation.size(); i++)
		{
			string option = options[i % options.size()];
			auto vertices = _triangulation[i].Vertices();
			if (vertices.size() == 3)
				script.PlotTriangle(vertices[0], vertices[1], vertices[2], option);
			else
				assert(false && "To implement");
		}
	}

	void ExportToMatlab(string color = "r") const override
	{
		MatlabScript script;
		script.PlotPolygon(_vertices, color);
		script.Out() << endl;
	}

	void ExportVerticesToMatlab(string color = "r") const
	{
		MatlabScript script;
		script.PlotPolygonEdges(_vertices, color);
		script.Out() << endl;
	}

	void ExportCGALPolyToMatlab(string color = "b") const
	{
		MatlabScript script;
		script.PlotPolygonEdges(CGALWrapper::ToVertices(_cgalPolygon), color);
		script.Out() << endl;
	}

	void ExportCGALPolyToMatlab(const CGAL::Polygon_2<chosenKernel>& cgalPoly, string color = "b") const
	{
		MatlabScript script;
		script.PlotPolygonEdges(CGALWrapper::ToVertices(cgalPoly), color);
		script.Out() << endl;
	}

	//------------------------------//
	//           Integral           //
	//------------------------------//

	vector<DomPoint> QuadraturePoints() const override
	{
		return _quadraturePoints;
	}

	double Integral(RefFunction boundingBoxDefinedFunction) const override
	{
		assert(_triangulation.size() > 0);

		// For the inner triangles, the bounding box is the domain
		DomFunction boundingBoxFunction = [this, boundingBoxDefinedFunction](const DomPoint& boundingBoxPoint)
		{
			RefPoint p = _boundingBox.ConvertToReference(boundingBoxPoint);
			return boundingBoxDefinedFunction(p);
		};

		double integral = 0;
		for (const Triangle& t : _triangulation)
			integral += t.Integral(boundingBoxFunction);
		return integral;
	}

	double Integral(RefFunction boundingBoxDefinedFunction, int polynomialDegree) const override
	{
		assert(_triangulation.size() > 0);

		// For the inner triangles, the bounding box is the domain
		DomFunction boundingBoxFunction = [this, boundingBoxDefinedFunction](const DomPoint& boundingBoxPoint)
		{
			RefPoint p = _boundingBox.ConvertToReference(boundingBoxPoint);
			return boundingBoxDefinedFunction(p);
		};

		double integral = 0;
		for (const Triangle& t : _triangulation)
			integral += t.Integral(boundingBoxFunction, polynomialDegree);
		return integral;
	}





	inline double DetJacobian(const RefPoint& pointInReferenceSquare) const
	{
		return _boundingBox.DetJacobian(pointInReferenceSquare);
	}
	inline DimMatrix<2> InverseJacobianTranspose(const RefPoint& pointInReferenceSquare) const
	{
		return _boundingBox.InverseJacobianTranspose(pointInReferenceSquare);
	}
	inline int DetJacobianDegree() const
	{
		return _boundingBox.DetJacobianDegree();
	}

	DomPoint ConvertToDomain(const RefPoint& pointInReferenceSquare) const
	{
		return _boundingBox.ConvertToDomain(pointInReferenceSquare);
	}

	RefPoint ConvertToReference(const DomPoint& domainPoint) const
	{
		return _boundingBox.ConvertToReference(domainPoint);
	}

	void Serialize(ostream& os) const override
	{
		os << "Polygon";
	}

	~Polygon()
	{}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	void UnitTests() const override
	{
		assert(_cgalPolygon.is_simple());
		assert(_cgalPolygon.is_counterclockwise_oriented());

		double sumMeasures = 0;
		for (const Triangle& t : _triangulation)
			sumMeasures += t.Measure();
		double eps = Utils::Eps*_measure;
		if (abs(sumMeasures - CGAL::to_double(_cgalPolygon.area())) >= eps)
		{
			this->ExportToMatlab();
			assert(abs(sumMeasures - CGAL::to_double(_cgalPolygon.area())) < eps);
		}

		PhysicalShape<2>::UnitTests();
	}

	static void Test()
	{
		DomPoint lowerLeft(-1, -1);
		DomPoint lowerRight(1, -1);
		DomPoint upperRight(1, 1);
		DomPoint upperLeft(-1, 1);
		vector<DomPoint> vertices{ lowerLeft, lowerRight, upperRight, upperLeft };
		Polygon polygRefSquare(vertices, true, false);

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
		RefFunction refX = [](const RefPoint& p) { return p.X; };
		double integral = polygRefSquare.Integral(refX);
		assert(integral == 0);

		//--------------------------------------//
		RefFunction anyFunction = [](const RefPoint& p) { return p.X + p.Y*p.Y + 1; };

		RectangleShape realRefSquare(&lowerLeft, 2);
		double integralOverRealSquare = realRefSquare.Integral(anyFunction);
		double integralOverPolygonalSquare = polygRefSquare.Integral(anyFunction);
		assert(abs(integralOverRealSquare - integralOverPolygonalSquare) < 1e-10);

		//--------------------------------------//
		double h = 2;
		lowerLeft = DomPoint(0, 0);
		lowerRight = DomPoint(h, 0);
		upperRight = DomPoint(h, h);
		upperLeft = DomPoint(0, h);
		Polygon polygSquare(vector<DomPoint>{ lowerLeft, lowerRight, upperRight, upperLeft }, true, false);
		RectangleShape realSquare(&lowerLeft, h);
		integralOverRealSquare = realSquare.Integral(anyFunction);
		integralOverPolygonalSquare = polygSquare.Integral(anyFunction);
		assert(abs(integralOverRealSquare - integralOverPolygonalSquare) < 1e-10);

		//--------------------------------------//
		h = 1;
		lowerLeft = DomPoint(0, 0);
		lowerRight = DomPoint(h, 0);
		upperRight = DomPoint(h, h);
		upperLeft = DomPoint(0, h);
		Polygon polygSquare1(vector<DomPoint>{ lowerLeft, lowerRight, upperRight, upperLeft }, true, false);
		realSquare = RectangleShape(&lowerLeft, h);
		integralOverRealSquare = realSquare.Integral(anyFunction);
		integralOverPolygonalSquare = polygSquare1.Integral(anyFunction);
		assert(abs(integralOverRealSquare - integralOverPolygonalSquare) < 1e-10);
	}
};