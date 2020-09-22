#pragma once
#include <cstdio>
#include<iostream>
using namespace std;

class Point
{
public:
	static double Tolerance;

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

double Point::Tolerance = 1e-10;

bool operator==(Point const& p1, Point const& p2)
{
	return abs(p1.X - p2.X) <= Point::Tolerance && abs(p1.Y - p2.Y) <= Point::Tolerance && abs(p1.Z - p2.Z) <= Point::Tolerance;
}
bool operator!=(Point const& p1, Point const& p2)
{
	return !(p1 == p2);
}
bool operator<(Point const& p1, Point const& p2)
{
	if (p1.X < p2.X)
		return true;
	if (p1.X == p2.X && p1.Y < p2.Y)
		return true;
	if (p1.X == p2.X && p1.Y == p2.Y && p1.Z < p2.Z)
		return true;
	return false;
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

	static DomPoint Middle(const DomPoint& A, const DomPoint& B)
	{
		return DomPoint((A.X + B.X) / 2, (A.Y + B.Y) / 2, (A.Z + B.Z) / 2);
	}
};