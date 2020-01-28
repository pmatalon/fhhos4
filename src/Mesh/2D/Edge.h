#pragma once
#include "../../DG/Poisson_DG_Face.h"
#include "../../HHO/Poisson_HHO_Face.h"
#include "../CartesianShape.h"
#include "../../Utils/Geometry.h"

class EdgeShape : public GeometricShapeWithConstantJacobian<1>
{
private:
	double _width;
	DomPoint _center;
public:
	Vertex* Vertex1;
	Vertex* Vertex2;

	EdgeShape(Vertex* v1, Vertex* v2) : GeometricShapeWithConstantJacobian<1>()
	{
		Vertex1 = v1;
		Vertex2 = v2;
		Init();
	}

	inline void Init()
	{
		_width = sqrt(pow(Vertex2->X - Vertex1->X, 2) + pow(Vertex2->Y - Vertex1->Y, 2));
		_center = DomPoint((Vertex1->X + Vertex2->X) / 2, (Vertex1->Y + Vertex2->Y) / 2);
	}

	inline ReferenceShape<1>* RefShape() const override
	{
		return &CartesianShape<2, 1>::RefCartShape;
	}

	inline vector<Vertex*> Vertices() const override
	{
		return vector<Vertex*> {Vertex1, Vertex2};
	}

	inline double Diameter() const override
	{
		return _width;
	}
	inline double Measure() const override
	{
		return _width;
	}
	inline DomPoint Center() const override
	{
		return _center;
	}
	inline bool Contains(DomPoint p) const override
	{
		return Geometry::IsInSegment(*Vertex1, *Vertex2, p);
	}

	inline double DetJacobian() const override
	{
		return _width / RefShape()->Measure();
	}
	inline DimMatrix<1> InverseJacobianTranspose() const override
	{
		assert(false);
	}

	DomPoint ConvertToDomain(RefPoint referenceElementPoint) const override
	{
		double x1 = Vertex1->X;
		double x2 = Vertex2->X;
		double y1 = Vertex1->Y;
		double y2 = Vertex2->Y;

		double t = referenceElementPoint.X;

		DomPoint p((x2 - x1) / 2 * t + (x2 + x1) / 2, (y2 - y1) / 2 * t + (y2 + y1) / 2);
		return p;
	}

	RefPoint ConvertToReference(DomPoint domainPoint) const override
	{
		double x1 = Vertex1->X;
		double x2 = Vertex2->X;
		double y1 = Vertex1->Y;
		double y2 = Vertex2->Y;

		double x = domainPoint.X;
		double y = domainPoint.Y;

		if (abs(x2 - x1) < abs(y2 - y1))
		{
			RefPoint p(2 / (y2 - y1) * y - (y2 + y1) / (y2 - y1));
			return p;
		}
		else
		{
			RefPoint p(2 / (x2 - x1) * x - (x2 + x1) / (x2 - x1));
			return p;
		}
	}

	void Serialize(ostream& os) const override
	{
		Vertex1->Serialize(os, 2);
		os << "--";
		Vertex2->Serialize(os, 2);
	}
};

class Edge : public Poisson_DG_Face<2>, public Poisson_HHO_Face<2>
{
private:
	EdgeShape* _shape;
public:

	Edge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1, Element<2>* element2) : 
		Poisson_DG_Face(number, element1, element2),
		Poisson_HHO_Face(number, element1, element2),
		Face(number, element1, element2)
	{
		_shape = new EdgeShape(v1, v2);
	}

	Edge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1) :
		Edge(number, v1, v2, element1, nullptr)
	{ }

	inline Vertex* Vertex1() const
	{
		return _shape->Vertex1;
	}
	inline Vertex* Vertex2() const
	{
		return _shape->Vertex2;
	}

	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	GeometricShapeWithReferenceShape<1>* Shape() const override
	{
		return _shape;
	}

	Face<2>* CreateSameGeometricFace(BigNumber number, Element<2>* element1)
	{
		Face<2>* copy = new Edge(number, _shape->Vertex1, _shape->Vertex2, element1);
		copy->IsDomainBoundary = this->IsDomainBoundary;
		return copy;
	}

	void ExportFaceToMatlab(FILE* file)
	{
		//             Number  x1    y1    x2    y2 IsDomainBoundary
		fprintf(file, "%lu %.17g %.17g %.17g %.17g %d\n", this->Number, _shape->Vertex1->X, _shape->Vertex1->Y, _shape->Vertex2->X, _shape->Vertex2->Y, this->IsDomainBoundary);
	}

	virtual ~Edge()
	{
		if (_shape)
			delete _shape;
	}
};