#pragma once
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
};

typedef Point RefPoint; // Reference point: X, Y, Z are in the reference interval [-1, 1]
typedef Point DomPoint; // Domain point