#pragma once
#include "Point.h"
#include "../Utils/Types.h"
using namespace std;

class Vertex : public DomPoint
{
public:
	BigNumber Number;

	Vertex(BigNumber number, double x) : DomPoint(x)
	{
		this->Number = number;
	}
	Vertex(BigNumber number, double x, double y) : DomPoint(x, y)
	{
		this->Number = number;
	}
	Vertex(BigNumber number, double x, double y, double z) : DomPoint(x, y, z)
	{
		this->Number = number;
	}
};

DimVector<2> operator-(Point const& A, Point const& B)
{
	DimVector<2> v;
	v << A.X - B.X, A.Y - B.Y;
	return v;
}

template<int Dim>
DimVector<Dim> Vect(Point A, Point B)
{
	DimVector<Dim> v;
	if (Dim == 1)
		v << B.X - A.X;
	else if (Dim == 2)
		v << B.X - A.X, B.Y - A.Y;
	else if (Dim == 3)
		v << B.X - A.X, B.Y - A.Y, B.Z - A.Z;
	return v;
}

template<int Dim>
DimVector<Dim> Vect(Vertex* A, Vertex* B)
{
	return Vect<Dim>(*A, *B);
}