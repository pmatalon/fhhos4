#pragma once
#include <functional>
#include <vector>
#include "GaussLegendre.h"

class Utils
{
public:
	static double Integral(int nPoints, std::function<double(double)> func, double x1, double x2)
	{
		GaussLegendre gs(nPoints);
		if (x1 == -1 && x2 == 1)
			return gs.Quadrature(func);
		else
			return gs.Quadrature(func, x1, x2);
	}
	
	static double Integral(std::function<double(double)> func, double x1, double x2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2);
	}

	static double Integral(int nPoints, std::function<double(double)> func, RefInterval interval)
	{
		return Integral(nPoints, func, interval.Left, interval.Right);
	}

	static double Integral(std::function<double(double)> func, RefInterval interval)
	{
		return Integral(func, interval.Left, interval.Right);
	}


	// Integral on [-1, 1] x [-1, 1]
	/*static double IntegralOnReferenceInterval(std::function<double(double, double)> func)
	{
		GaussLegendre gs;
		return gs.Quadrature(func);
	}*/

	// Integral on [x1, x2] x [y1, y2]
	static double Integral(int nPoints, std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		GaussLegendre gs(nPoints);
		if (x1 == -1 && x2 == 1 && y1 == -1 && y2 == 1)
			return gs.Quadrature(func);
		else
			return gs.Quadrature(func, x1, x2, y1, y2);
	}

	static double Integral(std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2);
	}

	static double Integral(int nPoints, std::function<double(double, double)> func, RefInterval xInterval, RefInterval yInterval)
	{
		return Integral(nPoints, func, xInterval.Left, xInterval.Right, yInterval.Left, yInterval.Right);
	}

	static double Integral(std::function<double(double, double)> func, RefInterval xInterval, RefInterval yInterval)
	{
		return Integral(func, xInterval.Left, xInterval.Right, yInterval.Left, yInterval.Right);
	}

	static int Binomial(int  n, int p)
	{
		if (p != 0 && n != p)
			return Binomial(n - 1, p) + Binomial(n - 1, p - 1);
		return 1;
	}

};