#pragma once
#include "GeometricShapeWithReferenceShape.h"
#include "2D/Quadrilateral.h"
using namespace std;

template <int Dim>
class AgglomerateShape : public GeometricShapeWithReferenceShape<Dim>
{
protected:
	vector<Vertex*> _vertices;

	double _diameter;
	double _measure;
	DomPoint _center;
	double _inRadius;

	vector<GeometricShapeWithReferenceShape<Dim>*> _subShapes;
	GeometricShapeWithReferenceShape<Dim>* _boundingBox = nullptr;

public:
	AgglomerateShape(GeometricShapeWithReferenceShape<Dim>* s) : GeometricShapeWithReferenceShape<Dim>()
	{
		Add(s);
	}

	void Add(GeometricShapeWithReferenceShape<Dim>* s)
	{
		AgglomerateShape* a = dynamic_cast<AgglomerateShape*>(s);
		if (a)
		{
			for (GeometricShapeWithReferenceShape<Dim>* ss : s->SubShapes())
			{
				GeometricShapeWithReferenceShape<Dim>* copy = ss->CreateCopy();
				assert(!this->Contains(copy->Center()));
				_subShapes.push_back(copy);
			}
		}
		else
		{
			GeometricShapeWithReferenceShape<Dim>* copy = s->CreateCopy();
			assert(!this->Contains(copy->Center()));
			_subShapes.push_back(copy);
		}
	}

	bool IsMadeOfSubShapes() const override
	{
		return true;
	}
	vector<GeometricShapeWithReferenceShape<Dim>*> SubShapes() const override
	{
		return _subShapes;
	}

	GeometricShapeWithReferenceShape<Dim>* CreateCopy() const
	{
		AgglomerateShape* copy = new AgglomerateShape(*this);
		copy->_boundingBox = this->_boundingBox->CreateCopy();
		return copy;
	}

	bool IsDegenerated() const override
	{
		return _subShapes.empty();
	}

	void ReshapeByMovingIntersection(Vertex* oldIntersect, Vertex* newIntersect) override
	{
		if (*oldIntersect == *newIntersect)
			return;

		vector<GeometricShapeWithReferenceShape<Dim>*> subShapesToDelete;
		for (GeometricShapeWithReferenceShape<Dim>* ss : _subShapes)
		{
			if (ss->HasVertex(oldIntersect))
			{
				if (ss->HasVertex(newIntersect)) // 2 vertices become one, the subshape loses a vertex
				{
					TriangleShape* triangle = dynamic_cast<TriangleShape*>(ss);
					if (triangle) // if it's a triangle, it degenerates
						subShapesToDelete.push_back(ss);
					else
					{
						QuadrilateralShape* quad = dynamic_cast<QuadrilateralShape*>(ss);
						if (quad) // the quadrilateral turns into a triangle
						{
							vector<Vertex*> vertices(3);
							for (Vertex* v : quad->Vertices())
							{
								if (v != oldIntersect)
									vertices.push_back(v);
							}
							ss = dynamic_cast<GeometricShapeWithReferenceShape<Dim>*>(new TriangleShape(vertices[0], vertices[1], vertices[2]));
							assert(ss);
						}
						else
						{
							AgglomerateShape<Dim>* agglo = dynamic_cast<AgglomerateShape<Dim>*>(ss);
							if (agglo)
							{
								agglo->ReshapeByMovingIntersection(oldIntersect, newIntersect);
								if (agglo->SubShapes().empty())
									subShapesToDelete.push_back(ss);
							}
							else
								assert(false);
						}
					}
				}
				else
				{
					CartesianShape<Dim>* cart = dynamic_cast<CartesianShape<Dim>*>(ss);
					if (cart)
						assert(false && "TODO create a quadrilateral");
					ss->ReshapeByMovingIntersection(oldIntersect, newIntersect);
					if (ss->IsDegenerated())
						subShapesToDelete.push_back(ss);
				}
			}
		}

		// Delete the shapes fully included in another one
		for (auto ss1 : _subShapes)
		{
			if (ss1->IsIn(subShapesToDelete))
				continue;
			for (auto ss2 : _subShapes)
			{
				if (ss2->IsIn(subShapesToDelete))
					continue;
				if (ss1 != ss2 && ss1->Contains(ss2))
					subShapesToDelete.push_back(ss2);
			}
		}

		// Actual deletion of the shapes that must be deleted
		for (auto sToDelete : subShapesToDelete)
		{
			int i = 0;
			for (i = 0; i < _subShapes.size(); i++)
			{
				if (_subShapes[i] == sToDelete)
					break;
			}
			assert(i < _subShapes.size());
			_subShapes.erase(_subShapes.begin() + i);
			delete sToDelete;
		}
	}

	virtual void Init(vector<Vertex*> vertices)
	{
		_vertices = vertices;
		if (_subShapes.size() == 1)
		{
			auto s = _subShapes[0];
			_diameter = s->Diameter();
			_measure = s->Measure();
			_center = s->Center();
			_inRadius = s->InRadius();
			_boundingBox = _subShapes[0];
		}
		else
		{
			_diameter = 0;
			double sumX = 0;
			double sumY = 0;
			double sumZ = 0;
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
				if (Dim == 3)
					sumZ += v1->Z;
			}

			_center = DomPoint(sumX / _vertices.size(), sumY / _vertices.size(), sumZ / _vertices.size());

			_measure = 0;
			for (GeometricShapeWithReferenceShape<Dim>* s : _subShapes)
				_measure += s->Measure();

			// TODO
			_inRadius = 0;

			if (Dim == 2)
				_boundingBox = dynamic_cast<GeometricShapeWithReferenceShape<Dim>*>(Geometry::CreateBoundingBox(_vertices));
			else
				assert(false && "To be implemented");
		}
	}

	virtual ReferenceShape<Dim>* RefShape() const override
	{
		assert(false && "Dim-specific function");
	}

	inline vector<Vertex*> Vertices() const override
	{
		return _vertices;
	}
	/*void SetVertices(vector<Vertex*> newVertices)
	{
		_vertices = newVertices;
	}*/

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
	bool IsConvex() const override
	{
		if (_vertices.size() < 4)
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
		return true;
	}
	inline double InRadius() const override
	{
		return _inRadius;
	}
	inline bool Contains(DomPoint p) const override
	{
		for (GeometricShapeWithReferenceShape<Dim>* s : _subShapes)
		{
			if (s->Contains(p))
				return true;
		}
		return false;
	}
	GeometricShapeWithReferenceShape<Dim>* ClosestSubShape(DomPoint p)
	{
		GeometricShapeWithReferenceShape<Dim>* closestSubShape = nullptr;
		double minDistance = 0;
		for (GeometricShapeWithReferenceShape<Dim>* s : _subShapes)
		{
			double distance = Vect<Dim>(s->Center(), p).norm();
			if (!closestSubShape || distance < minDistance)
			{
				closestSubShape = s;
				minDistance = distance;
			}
		}
		return closestSubShape;
	}

	void ExportSubShapesToMatlab() const override
	{
		MatlabScript script;
		vector<string> options = { "r", "b", "g", "c", "m", "y" };
		for (int i = 0; i < _subShapes.size(); i++)
		{
			string option = options[i % options.size()];
			/*RotatingList<Vertex*> vertices(_subShapes[i]->Vertices());
			for (int i = 0; i < vertices.Size(); i++)
			{
				auto v1 = vertices.GetAndMoveNext();
				auto v2 = vertices.Get();
				script.PlotSegment(v1, v2, option);
			}*/
			auto vertices = _subShapes[i]->Vertices();
			if (vertices.size() == 3)
				script.PlotTriangle(*vertices[0], *vertices[1], *vertices[2], option);
			else
				assert(false && "To implement");
		}
	}

	void CheckIfTwoSubShapesOverlap() const
	{
		for (auto ss1 : _subShapes)
		{
			for (auto ss2 : _subShapes)
			{
				if (ss1 == ss2)
					continue;

				bool overlap = ss1->Contains(ss2->Center());
				if (overlap)
				{
					ExportSubShapesToMatlab();
					bool b = ss1->Contains(ss2);
				}
				assert(!overlap && "Two subshapes overlap");
			}
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

		QuadrilateralShape* quad = dynamic_cast<QuadrilateralShape*>(_boundingBox); // TODO remove

		double integral = 0;
		for (GeometricShapeWithReferenceShape<Dim>* t : _subShapes)
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
		for (GeometricShapeWithReferenceShape<Dim>* s : _subShapes)
			integral += s->Integral(boundingBoxFunction, polynomialDegree);
		return integral;
	}

	inline double DetJacobian(RefPoint pointInReferenceSquare) const
	{
		return _boundingBox->DetJacobian(pointInReferenceSquare);
	}
	inline DimMatrix<Dim> InverseJacobianTranspose(RefPoint pointInReferenceSquare) const
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
		os << "Agglomerate";
	}

	virtual ~AgglomerateShape()
	{
		for (GeometricShapeWithReferenceShape<Dim>* s : _subShapes)
			delete s;
		if (_subShapes.size() > 1)
			delete _boundingBox;
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	void UnitTests() const
	{
		GeometricShapeWithReferenceShape<Dim>::UnitTests();

		assert(_subShapes.size() > 0);
		CheckIfTwoSubShapesOverlap();
	}
};

template<>
ReferenceShape<2>* AgglomerateShape<2>::RefShape() const
{
	if (_subShapes.size() == 1)
		return _subShapes[0]->RefShape();
	else
		return &RectangleShape::RefCartShape;
}