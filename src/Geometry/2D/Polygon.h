#pragma once
#include "../../Utils/RotatingList.h"
#include "../../Utils/MatlabScript.h"
#include "../CartesianShape.h"
#include "../CGALConvert.h"
#include "Triangle.h"
#include "Quadrilateral.h"
#include <CGAL/convex_hull_2.h>
#include <cassert>
#include <list>
using namespace std;

typedef exactKernel chosenKernel;

class Polygon : public PhysicalShape<2>
{
private:
	vector<Vertex*> _vertices;
	//vector<Vertex*> _nonColinearVertices;

	CGAL::Polygon_2<chosenKernel> _cgalPolygon; // Polygon of the CGAL library

	double _diameter;
	double _measure;
	Vertex* _center;
	double _inRadius;

	vector<PhysicalShape<2>*> _triangulation;
	Quadrilateral* _boundingBox = nullptr;
	vector<DomPoint> _quadraturePoints;

public:
	Polygon(const vector<Vertex*>& vertices, bool createTriangulationAndBoundingBox, bool checkForColinearVertices)
		: _vertices(vertices)
	{
		assert(vertices.size() >= 3);
		InitCGALPolygon(checkForColinearVertices);
		Init(createTriangulationAndBoundingBox);
	}

	Polygon(const Polygon& shape) = default;

private:
	Polygon(const CGAL::Polygon_2<chosenKernel>& cgalPolygon, bool createTriangulationAndBoundingBox)
		: _cgalPolygon(cgalPolygon)
	{
		_vertices = CGALWrapper::ToVertices(cgalPolygon);
		Init(createTriangulationAndBoundingBox);
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
			ExportToMatlab("m");
			Utils::FatalError("Non simple polygon detected.");
		}
		else if (!_cgalPolygon.is_counterclockwise_oriented())
		{
			ExportToMatlab("m");
			Utils::FatalError("Clockwise-oriented polygon detected.");
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
			_triangulation = CGALTriangulation();

		assert(_triangulation.size() > 0);

		//InitQuadraturePoints();
	}

	void ComputeBoundingBox()
	{
		if (_boundingBox)
			return;
		_boundingBox = CreateBoundingBox(_vertices);
	}

	static Quadrilateral* CreateBoundingBox(const vector<Vertex*>& vertices)
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

		Quadrilateral* boundingRectangle = new Quadrilateral(lowerLeft, lowerRight, upperRight, upperLeft);
		return boundingRectangle;
	}

private:
	static vector<PhysicalShape<2>*> ConvexTriangulation(const vector<Vertex*>& vertices)
	{
		assert(vertices.size() > 2);
		vector<PhysicalShape<2>*> triangles;
		if (vertices.size() == 3)
		{
			Triangle* triangle = new Triangle(vertices[0], vertices[1], vertices[2]);
			triangles.push_back(triangle);
		}
		else if (vertices.size() == 4)
		{
			Triangle* triangle1 = new Triangle(vertices[0], vertices[1], vertices[2]);
			triangles.push_back(triangle1);
			Triangle* triangle2 = new Triangle(vertices[2], vertices[3], vertices[0]);
			triangles.push_back(triangle2);
			// Uncomment when Contains() is implemented for the Quadrilateral
			//Quadrilateral* q = new Quadrilateral(vertices[0], vertices[1], vertices[2], vertices[3]);
			//triangles.push_back(q);
		}
		else
		{
			Triangle* triangle = new Triangle(vertices[0], vertices[1], vertices[2]);
			triangles.push_back(triangle);

			vector<Vertex*> remainingVertices;
			for (int i = 2; i < vertices.size(); i++)
				remainingVertices.push_back(vertices[i]);
			remainingVertices.push_back(vertices[0]);
			vector<PhysicalShape<2>*> otherShapes = ConvexTriangulation(remainingVertices);
			for (auto s : otherShapes)
				triangles.push_back(s);
		}
		return triangles;
	}
	static vector<PhysicalShape<2>*> BarycentricTriangulation(const vector<Vertex*>& vertices)
	{
		// Requirement: the polygon defined by the vertices must be convex!

		vector<PhysicalShape<2>*> triangles;

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

	/*static void my_cgal_failure_handler(
		const char *type,
		const char *expr,
		const char* file,
		int line,
		const char* msg)
	{
		// report the error in some way.
	}*/

	vector<PhysicalShape<2>*> CGALTriangulation()
	{
		vector<PhysicalShape<2>*> triangulation;

		list<CGAL::Partition_traits_2<chosenKernel>::Polygon_2> partition_polys;

		//CGAL::Failure_function prev; // save the CGAL failure function
		//prev = CGAL::set_error_handler(my_cgal_failure_handler); // replace it with my own
		CGAL::greene_approx_convex_partition_2(_cgalPolygon.vertices_begin(), _cgalPolygon.vertices_end(), back_inserter(partition_polys));
		//CGAL::set_error_handler(prev); // put the old one back

		for (CGAL::Partition_traits_2<chosenKernel>::Polygon_2 p : partition_polys)
		{
			vector<Vertex*> vertices = CGALWrapper::ToVertices(p);
			vector<PhysicalShape<2>*> subTriangles = ConvexTriangulation(vertices);
			for (auto tri : subTriangles)
				triangulation.push_back(tri);
		}

		return triangulation;
	}

	void InitQuadraturePoints()
	{
		for (PhysicalShape<2>* t : _triangulation)
		{
			for (const RefPoint& refPoint : t->RefShape()->QuadraturePoints())
			{
				DomPoint domPoint = t->ConvertToDomain(refPoint);//this->ConvertToDomainAndSaveResult(refPoint);
				_quadraturePoints.push_back(domPoint);
			}
		}
	}

public:
	bool ConvexHullEmbeds(PhysicalShape<2>* s) const override
	{
		vector<CGAL::Point_2<chosenKernel>> convexHullPoints;
		CGAL::convex_hull_2(_cgalPolygon.vertices_begin(), _cgalPolygon.vertices_end(), back_inserter(convexHullPoints));

		CGAL::Polygon_2<chosenKernel> convexHull(convexHullPoints.begin(), convexHullPoints.end());

		for (Vertex* v : s->Vertices())
		{
			if (!CGALWrapper::CGALPolyContains(convexHull, CGAL::Point_2<chosenKernel>(v->X, v->Y)))
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

	vector<PhysicalShape<2>*> IntersectionWith(PhysicalShape<2>* other) override
	{
		// Compute the intersection
		vector<PhysicalShape<2>*> intersection;
		list<CGAL::Polygon_with_holes_2<exactKernel>> intersectionPolygons;
		Polygon* otherPolygon = dynamic_cast<Polygon*>(other);
		if (otherPolygon)
			CGAL::intersection(_cgalPolygon, otherPolygon->_cgalPolygon, back_inserter(intersectionPolygons));
		else
		{
			auto otherShape = CGALWrapper::CreatePolygon<exactKernel>(other->Vertices());
			CGAL::intersection(_cgalPolygon, otherShape, back_inserter(intersectionPolygons));
		}

		// The intersection can be made of multiple polygons
		for (auto cgalPolyWithHoles : intersectionPolygons)
		{
			CGAL::Polygon_2<chosenKernel> cgalPoly = cgalPolyWithHoles.outer_boundary();

			if (cgalPoly.area() < Utils::NumericalZero * this->Measure())
				continue; // the intersection is the interface

			/*if (!cgalPoly.is_simple())
			{
				cout << "% Poly 1" << endl;
				this->ExportToMatlab("r");
				cout << "% Poly 2" << endl;
				other->ExportToMatlab("b");
				cout << "% Intersection" << endl;
				ExportCGALPolyToMatlab(cgalPoly, "m");
				Utils::Warning("Intersection polygon is not simple.");
			}*/
			Polygon* intersectionPolygon = new Polygon(cgalPoly, true);
			intersection.push_back(intersectionPolygon);
		}
		return intersection;
	}


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
	vector<PhysicalShape<2>*> SubShapes() const override
	{
		assert(_triangulation.size() > 0);
		return _triangulation;
	}

	virtual ReferenceShape<2>* RefShape() const override
	{
		return &RectangleShape::RefCartShape;
	}

	inline const vector<Vertex*>& Vertices() const override
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
			auto vertices = _triangulation[i]->Vertices();
			if (vertices.size() == 3)
				script.PlotTriangle(*vertices[0], *vertices[1], *vertices[2], option);
			else
				assert(false && "To implement");
		}
	}

	void ExportToMatlab(string color = "r") const override
	{
		//cout << "%-- Shape using vertices:" << endl;
		ExportVerticesToMatlab(color);
		/*cout << "%-- Shape using CGAL vertices:" << endl;
		ExportCGALPolyToMatlab();
		if (_triangulation.size() > 0)
		{
			cout << "%-- Subshapes:" << endl;
			ExportSubShapesToMatlab();
		}*/
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
			RefPoint p = _boundingBox->ConvertToReference(boundingBoxPoint);
			return boundingBoxDefinedFunction(p);
		};

		double integral = 0;
		for (PhysicalShape<2>* t : _triangulation)
			integral += t->Integral(boundingBoxFunction);
		return integral;
	}

	double Integral(RefFunction boundingBoxDefinedFunction, int polynomialDegree) const override
	{
		assert(_triangulation.size() > 0);

		// For the inner triangles, the bounding box is the domain
		DomFunction boundingBoxFunction = [this, boundingBoxDefinedFunction](const DomPoint& boundingBoxPoint)
		{
			RefPoint p = _boundingBox->ConvertToReference(boundingBoxPoint);
			return boundingBoxDefinedFunction(p);
		};

		double integral = 0;
		for (PhysicalShape<2>* t : _triangulation)
			integral += t->Integral(boundingBoxFunction, polynomialDegree);
		return integral;
	}





	inline double DetJacobian(const RefPoint& pointInReferenceSquare) const
	{
		return _boundingBox->DetJacobian(pointInReferenceSquare);
	}
	inline DimMatrix<2> InverseJacobianTranspose(const RefPoint& pointInReferenceSquare) const
	{
		return _boundingBox->InverseJacobianTranspose(pointInReferenceSquare);
	}
	inline int DetJacobianDegree() const
	{
		return _boundingBox->DetJacobianDegree();
	}

	DomPoint ConvertToDomain(const RefPoint& pointInReferenceSquare) const
	{
		return _boundingBox->ConvertToDomain(pointInReferenceSquare);
	}

	RefPoint ConvertToReference(const DomPoint& domainPoint) const
	{
		return _boundingBox->ConvertToReference(domainPoint);
	}

	void Serialize(ostream& os) const override
	{
		os << "Polygon";
	}

	~Polygon()
	{
		for (PhysicalShape<2>* t : _triangulation)
			delete t;
		delete _center;
		delete _boundingBox;
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	void UnitTests() const override
	{
		assert(_cgalPolygon.is_simple());
		assert(_cgalPolygon.is_counterclockwise_oriented());

		double sumMeasures = 0;
		for (PhysicalShape<2>* t : _triangulation)
			sumMeasures += t->Measure();
		double eps = 1e-4*_measure;
		if (abs(sumMeasures - _cgalPolygon.area()) >= eps)
		{
			this->ExportToMatlab();
			assert(abs(sumMeasures - _cgalPolygon.area()) < eps);
		}

		PhysicalShape<2>::UnitTests();
	}

	static void Test()
	{
		int number = 0;
		Vertex lowerLeft(number, -1, -1);
		Vertex lowerRight(number, 1, -1);
		Vertex upperRight(number, 1, 1);
		Vertex upperLeft(number, -1, 1);
		vector<Vertex*> vertices{ &lowerLeft, &lowerRight, &upperRight, &upperLeft };
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
		lowerLeft = Vertex(number, 0, 0);
		lowerRight = Vertex(number, h, 0);
		upperRight = Vertex(number, h, h);
		upperLeft = Vertex(number, 0, h);
		Polygon polygSquare(vector<Vertex*>{ &lowerLeft, &lowerRight, &upperRight, &upperLeft }, true, false);
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
		Polygon polygSquare1(vector<Vertex*>{ &lowerLeft, &lowerRight, &upperRight, &upperLeft }, true, false);
		realSquare = RectangleShape(&lowerLeft, h);
		integralOverRealSquare = realSquare.Integral(anyFunction);
		integralOverPolygonalSquare = polygSquare1.Integral(anyFunction);
		assert(abs(integralOverRealSquare - integralOverPolygonalSquare) < 1e-10);
	}
};