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

typedef Point RefPoint; // Reference point: X, Y, Z are in the reference interval [-1, 1]
typedef Point DomPoint; // Domain point