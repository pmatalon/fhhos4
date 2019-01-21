#pragma once
#include <functional>
#include <vector>
#include "GaussLegendre.h"

typedef unsigned int BigNumber;

class Utils
{
public:
	//-------------//
	// Integral 1D //
	//-------------//

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

	/*static double Integral(int nPoints, std::function<double(double)> func, DefInterval interval)
	{
		return Integral(nPoints, func, interval.Left, interval.Right);
	}

	static double Integral(std::function<double(double)> func, DefInterval interval)
	{
		return Integral(func, interval.Left, interval.Right);
	}*/

	//-------------//
	// Integral 2D //
	//-------------//

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

	/*static double Integral(int nPoints, std::function<double(double, double)> func, DefInterval xInterval, DefInterval yInterval)
	{
		return Integral(nPoints, func, xInterval.Left, xInterval.Right, yInterval.Left, yInterval.Right);
	}

	static double Integral(std::function<double(double, double)> func, DefInterval xInterval, DefInterval yInterval)
	{
		return Integral(func, xInterval.Left, xInterval.Right, yInterval.Left, yInterval.Right);
	}*/

	//-------------//
	// Integral 3D //
	//-------------//

	// Integral on [x1, x2] x [y1, y2]
	static double Integral(int nPoints, std::function<double(double, double, double)> func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		GaussLegendre gs(nPoints);
		if (x1 == -1 && x2 == 1 && y1 == -1 && y2 == 1 && z1 == -1 && z2 == 1)
			return gs.Quadrature(func);
		else
			return gs.Quadrature(func, x1, x2, y1, y2, z1, z2);
	}

	static double Integral(std::function<double(double, double, double)> func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2, z1, z2);
	}

	/*static double Integral(int nPoints, std::function<double(double, double, double)> func, DefInterval xInterval, DefInterval yInterval, DefInterval zInterval)
	{
		return Integral(nPoints, func, xInterval.Left, xInterval.Right, yInterval.Left, yInterval.Right, zInterval.Left, zInterval.Right);
	}

	static double Integral(std::function<double(double, double, double)> func, DefInterval xInterval, DefInterval yInterval, DefInterval zInterval)
	{
		return Integral(func, xInterval.Left, xInterval.Right, yInterval.Left, yInterval.Right, zInterval.Left, zInterval.Right);
	}*/

	//----------------------//
	// Binomial coefficient //
	//----------------------//

	static int Binomial(int n, int p)
	{
		if (p != 0 && n != p)
			return Binomial(n - 1, p) + Binomial(n - 1, p - 1);
		return 1;
	}

	static int Factorial(int n)
	{
		int f = 1;
		for (int i = 2; i <= n; i++)
			f *= i;
		return f;
	}

};