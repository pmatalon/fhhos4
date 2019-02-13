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