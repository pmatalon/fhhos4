#pragma once
#include "../Geometry/Point.h"
#include "../Utils/Utils.h"
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
	Vertex(BigNumber number, const DomPoint& p) : DomPoint(p)
	{
		this->Number = number;
	}

	bool IsIn(const vector<Vertex*>& list)
	{
		for (Vertex* v : list)
		{
			if (v == this)
				return true;
		}
		return false;
	}

	Vertex* GetClosestIn(const vector<Vertex*>& list)
	{
		double minDistance = 0;
		Vertex* closest = nullptr;
		for (auto v : list)
		{
			double distance = sqrt(pow(this->X - v->X, 2) + pow(this->Y - v->Y, 2) + pow(this->Z - v->Z, 2));
			if (!closest || distance < minDistance)
			{
				closest = v;
				minDistance = distance;
			}
		}
		return closest;
	}

	static vector<DomPoint> ToDomPoints(const vector<Vertex*>& list)
	{
		vector<DomPoint> points;
		for (auto v : list)
			points.push_back(*v);
		return points;
	}

	virtual ~Vertex() {}
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

template<int Dim>
DimVector<Dim> Vect(Vertex* A, Point B)
{
	return Vect<Dim>(*A, B);
}

template<int Dim>
DomPoint Middle(const DomPoint& A, const DomPoint& B)
{
	if (Dim == 1)
		return DomPoint((A.X + B.X) / 2);
	else if (Dim == 2)
		return DomPoint((A.X + B.X) / 2, (A.Y + B.Y) / 2);
	else if (Dim == 3)
		return DomPoint((A.X + B.X) / 2, (A.Y + B.Y) / 2, (A.Z + B.Z) / 2);
	assert(false);
}

template<int Dim>
DomPoint Middle(Vertex* A, Vertex* B)
{
	return Middle<Dim>(*A, *B);
}

template<int Dim>
bool AreColinear(const DomPoint& p1, const DomPoint& p2, const DomPoint& p3)
{
	DimVector<Dim> v1 = Vect<Dim>(p1, p2);
	DimVector<Dim> v2 = Vect<Dim>(p2, p3);
	return AreCollinear(v1, v2);
}

template<int Dim>
bool AreCollinear(const DimVector<Dim>& v1, const DimVector<Dim>& v2)
{
	double v1_norm_v2_norm = v1.norm()*v2.norm();
	return abs(abs(v1.dot(v2)) - v1_norm_v2_norm) / v1_norm_v2_norm < Utils::NumericalZero;
}

template<int Dim>
DomPoint AddVector(const DomPoint& p, const DimVector<Dim>& v)
{
	if (Dim == 1)
		return DomPoint(p.X + v[0]);
	else if (Dim == 2)
		return DomPoint(p.X + v[0], p.Y + v[1]);
	else if (Dim == 3)
		return DomPoint(p.X + v[0], p.Y + v[1], p.Z + v[2]);
}