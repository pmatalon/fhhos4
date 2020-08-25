#pragma once
#include "../../Mesh/Vertex.h"
#include "../PhysicalShapeWithConstantJacobian.h"
#include "../CartesianShape.h"
using namespace std;

class Segment : public PhysicalShapeWithConstantJacobian<1>
{
private:
	double _width;
	DomPoint _center;
public:
	Vertex* Vertex1;
	Vertex* Vertex2;

	Segment(Vertex* v1, Vertex* v2) : PhysicalShapeWithConstantJacobian<1>()
	{
		Vertex1 = v1;
		Vertex2 = v2;
		Init();
	}

	Segment(const Segment& shape) = default;

	inline void Init()
	{
		_width = sqrt(pow(Vertex2->X - Vertex1->X, 2) + pow(Vertex2->Y - Vertex1->Y, 2));
		_center = DomPoint((Vertex1->X + Vertex2->X) / 2, (Vertex1->Y + Vertex2->Y) / 2);
	}

	PhysicalShape<1>* CreateCopy() const
	{
		return new Segment(*this);
	}

	inline ReferenceShape<1>* RefShape() const override
	{
		return &CartesianShape<2, 1>::RefCartShape;
	}

	inline vector<Vertex*> Vertices() const override
	{
		return vector<Vertex*> {Vertex1, Vertex2};
	}

	bool IsDegenerated() const override
	{
		return *Vertex1 == *Vertex2;
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
	inline bool IsConvex() const override
	{
		return true;
	}
	inline double InRadius() const override
	{
		return 0;
	}
	inline bool Contains(DomPoint p) const override
	{
		return SegmentContains(*Vertex1, *Vertex2, p);
	}

	static bool SegmentContains(DomPoint A, DomPoint B, DomPoint P)
	{
		DimVector<2> AB = B - A;
		DimVector<2> AP = P - A;
		double AB_dot_AP = AB.dot(AP);
		return AB_dot_AP > 0 && AB_dot_AP < AB.dot(AB) && abs(AB_dot_AP - AB.norm()*AP.norm()) < Point::Tolerance;
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

	void ExportToMatlab(string color = "r") const override
	{
		MatlabScript script;
		script.PlotSegment(Vertex1, Vertex2, color);
	}

	void Serialize(ostream& os) const override
	{
		Vertex1->Serialize(os, 2);
		os << "--";
		Vertex2->Serialize(os, 2);
	}
};