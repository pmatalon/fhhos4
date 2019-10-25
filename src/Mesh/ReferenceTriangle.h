#pragma once
#include "ReferenceElement.h"
#include "../Utils/Dunavant/Dunavant.h"

/*struct BarycentricPoint
{
	double Alpha = 0;
	double Beta = 0;
	double Gamma = 0;
};

class TriangleStruct
{
private:
	Point A;
	Point B;
	Point C;
	double det;
public:
	TriangleStruct(Point a, Point b, Point c) : A(a), B(b), C(c)
	{
		det = (B.Y - C.Y)*(A.X - C.X) + (C.X - B.X)*(A.Y - C.Y);
	}

	BarycentricPoint ToBarycentric(Point p) const
	{
		BarycentricPoint bary;
		bary.Alpha = ((B.Y - C.Y)*(p.X - C.X) + (C.X - B.X)*(p.Y - C.Y)) / det;
		bary.Beta = ((C.Y - A.Y)*(p.X - C.X) + (A.X - C.X)*(p.Y - C.Y)) / det;
		bary.Gamma = 1 - bary.Alpha - bary.Beta;
		return bary;
	}

	Point ToCartesian(BarycentricPoint bary) const
	{
		Point p;
		p.X = bary.Alpha * A.X + bary.Beta * B.X + bary.Gamma * C.X;
		p.Y = bary.Alpha * A.Y + bary.Beta * B.Y + bary.Gamma * C.Y;
	}
};*/

class ReferenceTriangle : public ReferenceElement<2>
{
private:
	RefPoint A;
	RefPoint B;
	RefPoint C;
	//TriangleStruct _baryConverter;
	double _measure;

public:
	ReferenceTriangle() : 
		ReferenceElement<2>(),
		A(0, 0), B(1, 0), C(0, 1)//, _baryConverter(A, B, C)
	{
		_measure = 0.5 * abs(A.X * (B.Y - C.Y) + B.X * (C.Y - A.Y) + C.X * (A.Y - B.Y));
	}

	inline double Measure()
	{
		return _measure;
	}

	/*inline BarycentricPoint ConvertToBarycentric(RefPoint p) const
	{
		return _baryConverter.ToBarycentric(p);
	}
	inline RefPoint ConvertToRef(BarycentricPoint bary) const
	{
		return _baryConverter.ToCartesian(bary);
	}*/

	double ComputeIntegral(RefFunction func) const override
	{
		//return _measure * func(RefPoint(0, 0)); // TODO remove this!!!
		// TODO: change 1000
		Dunavant dunavant(1000);
		return _measure * dunavant.Quadrature(func);
	}
	double ComputeIntegral(RefFunction func, int polynomialDegree) const override
	{
		if (polynomialDegree == 0)
			return _measure * func(RefPoint(0, 0));
		Dunavant dunavant(polynomialDegree);
		return _measure * dunavant.Quadrature(func);
	}
	double ComputeIntegral(BasisFunction<2>* phi) const
	{
		return ReferenceElement<2>::ComputeIntegral(phi);
	}

};