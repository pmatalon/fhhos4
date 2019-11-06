#pragma once
#include "Point.h"
#include "../Utils/Types.h"

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

DimVector<2> Vect(Vertex* A, Vertex* B)
{
	DimVector<2> v;
	v << B->X - A->X, B->Y - A->Y;
	return v;
}