#pragma once
using namespace std;

class Point
{
public:
	double X = 0;
	double Y = 0;
	double Z = 0;

	Point() {}
	Point(double x)
	{
		this->X = x;
	}
	Point(double x, double y)
	{
		this->X = x;
		this->Y = y;
	}
	Point(double x, double y, double z)
	{
		this->X = x;
		this->Y = y;
		this->Z = z;
	}

	void Serialize(ostream& os, int dim) const
	{
		os << "(";
		if (dim == 1)
			os << X;
		else if (dim == 2)
			os << X << ", " << Y;
		else if (dim == 3)
			os << X << ", " << Y << ", " << Z;
		os << ")";
	}
};

bool operator==(Point const& p1, Point const& p2)
{
	return p1.X == p2.X && p1.Y == p2.Y && p1.Z == p2.Z;
}
bool operator!=(Point const& p1, Point const& p2)
{
	return p1.X != p2.X || p1.Y != p2.Y || p1.Z != p2.Z;
}


// Reference point: X, Y, Z are in the reference interval [-1, 1]
class RefPoint : public Point 
{
public:
	RefPoint() : Point() {}
	RefPoint(double x) : Point(x) {}
	RefPoint(double x, double y) : Point(x, y) {}
	RefPoint(double x, double y, double z) : Point(x, y, z) {}
};

// Domain point
class DomPoint : public Point
{
public:
	DomPoint() : Point() {}
	DomPoint(double x) : Point(x) {}
	DomPoint(double x, double y) : Point(x, y) {}
	DomPoint(double x, double y, double z) : Point(x, y, z) {}
};