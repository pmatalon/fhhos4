#pragma once
#include "../../Utils/Geometry.h"
#include "../CartesianShape.h"
using namespace std;

class PolygonalShape : public GeometricShapeWithReferenceShape<2>
{
private:
	vector<Vertex*> _vertices;

	double _diameter;
	double _measure;
	Vertex* _center;
	double _inRadius;

	vector<TriangleShape*> _triangulation;
	QuadrilateralShape* _boundingBox;

public:

	PolygonalShape(vector<Vertex*> vertices) 
		: _vertices(vertices)
	{
		assert(vertices.size() >= 3);
		Init();
	}

	void Init()
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

		_triangulation = Geometry::Triangulation(_vertices);

		_boundingBox = Geometry::CreateBoundingBox(_vertices);

		_measure = 0;
		for (TriangleShape* t : _triangulation)
			_measure += t->Measure();

		_inRadius = 0;
	}

	ReferenceShape<2>* RefShape() const
	{
		return &RectangleShape::RefCartShape;
	}

	inline vector<Vertex*> Vertices() const override
	{
		return _vertices;
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
	inline double InRadius() const override
	{
		return _inRadius;
	}
	inline bool Contains(DomPoint p) const override
	{
		assert(false && "Not implemented");
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
		for (TriangleShape* t : _triangulation)
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
		for (TriangleShape* t : _triangulation)
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
		for (TriangleShape* t : _triangulation)
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