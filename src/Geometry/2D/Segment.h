#pragma once
#include "../../Mesh/Vertex.h"
#include "../PhysicalShapeWithConstantJacobian.h"
#include "../CartesianShape.h"
using namespace std;

class Segment : public PhysicalShapeWithConstantJacobian<1>
{
private:
	DomPoint v1;
	DomPoint v2;
	double _width;
	DomPoint _center;
public:
	Segment() {}

	Segment(const DomPoint& p1, const DomPoint& p2)
		: PhysicalShapeWithConstantJacobian<1>(), v1(p1), v2(p2)
	{
		Init();
	}

	Segment(const Segment& shape) = default;

	inline DomPoint Vertex1() const { return v1; }
	inline DomPoint Vertex2() const { return v2; }

	inline void Init()
	{
		_width = sqrt(pow(v2.X - v1.X, 2) + pow(v2.Y - v1.Y, 2));
		_center = DomPoint((v1.X + v2.X) / 2, (v1.Y + v2.Y) / 2);
	}

	PhysicalShape<1>* CreateCopy() const
	{
		return new Segment(*this);
	}

	inline ReferenceShape<1>* RefShape() const override
	{
		return &CartesianShape<2, 1>::RefCartShape;
	}

	inline vector<DomPoint> Vertices() const override
	{
		return vector<DomPoint>{ v1, v2 };
	}

	bool IsDegenerated() const override
	{
		return v1 == v2;
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
	inline DomPoint InteriorPoint() const override
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
	inline bool Contains(const DomPoint& p) const override
	{
		return SegmentContains(v1, v2, p);
	}

	static bool SegmentContains(const DomPoint& A, const DomPoint& B, const DomPoint& P)
	{
		DimVector<2> AB = B - A;
		DimVector<2> AP = P - A;
		double AB_dot_AP = AB.dot(AP);
		return AB_dot_AP > 0 && AB_dot_AP < AB.dot(AB) && abs(AB_dot_AP - AB.norm()*AP.norm()) < Utils::Eps*AB_dot_AP;
	}

	inline double DetJacobian() const override
	{
		return _width / RefShape()->Measure();
	}
	inline DimMatrix<1> InverseJacobianTranspose() const override
	{
		assert(false);
	}

	DomPoint ConvertToDomain(const RefPoint& referenceElementPoint) const override
	{
		double x1 = v1.X;
		double x2 = v2.X;
		double y1 = v1.Y;
		double y2 = v2.Y;

		double t = referenceElementPoint.X;

		DomPoint p((x2 - x1) / 2 * t + (x2 + x1) / 2, (y2 - y1) / 2 * t + (y2 + y1) / 2);
		return p;
	}

	RefPoint ConvertToReference(const DomPoint& domainPoint) const override
	{
		double x1 = v1.X;
		double x2 = v2.X;
		double y1 = v1.Y;
		double y2 = v2.Y;

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
		script.PlotSegment(v1, v2, color);
	}

	void Serialize(ostream& os) const override
	{
		v1.Serialize(os, 2);
		os << "--";
		v2.Serialize(os, 2);
	}
};
