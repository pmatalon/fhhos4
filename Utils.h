#pragma once
#include <functional>
#include <vector>
#include "GaussLegendre.h"

class Utils
{
public:
	static double Integral(std::function<double(double)> func, double x1, double x2)
	{
		GaussLegendre gs;
		return gs.Quadrature(func, x1, x2);
	}
	static double Integral(int nPoints, std::function<double(double)> func, double x1, double x2)
	{
		GaussLegendre gs(nPoints);
		return gs.Quadrature(func, x1, x2);
	}

	// Integral on [-1, 1] x [-1, 1]
	static double IntegralOnReferenceInterval(std::function<double(double, double)> func)
	{
		GaussLegendre gs;
		return gs.Quadrature(func);
	}

	// Integral on [x1, x2] x [y1, y2]
	static double Integral(std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		GaussLegendre* gs = new GaussLegendre();
		if (x1 == x2)
		{
			std::function<double(double)> func1D = [func, x1](double y) {
				return func(x1, y);
			};
			return gs->Quadrature(func1D, y1, y2);
		}
		if (y1 == y2)
		{
			std::function<double(double)> func1D = [func, y1](double x) {
				return func(x, y1);
			};
			return gs->Quadrature(func1D, x1, x2);
		}
		return gs->Quadrature(func, x1, x2, y1, y2);
	}

	static int Binomial(int  n, int p)
	{
		if (p != 0 && n != p)
			return Binomial(n - 1, p) + Binomial(n - 1, p - 1);
		return 1;
	}

};