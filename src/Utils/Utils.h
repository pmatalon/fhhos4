#pragma once
#include <functional>
#include <vector>
#include "GaussLegendre.h"
#include "../Mesh/Point.h"
#include "../FunctionalBasis/BasisFunction.h"

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

	static double Integral(std::function<double(Point)> func, double x1, double x2)
	{
		function<double(double)> funcToIntegrate = [func](double x) {
			return func(x);
		};
		return Integral(funcToIntegrate, x1, x2);
	}

	static double Integral(BasisFunction<1>* phi, double x1, double x2)
	{
		function<double(double)> func = [phi](double x) {
			return phi->Eval(x);
		};
		int nPoints = (int)ceil(((double)phi->GetDegree() + 1)/2);
		return Utils::Integral(nPoints, func, x1, x2);
	}

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

	static double Integral(std::function<double(Point)> func, double x1, double x2, double y1, double y2)
	{
		function<double(double, double)> funcToIntegrate = [func](double x, double y) {
			Point p(x, y);
			return func(p);
		};
		return Integral(funcToIntegrate, x1, x2, y1, y2);
	}

	static double Integral(BasisFunction<2>* phi, double x1, double x2, double y1, double y2)
	{
		function<double(double, double)> func = [phi](double x, double y) {
			return phi->Eval(Point(x, y));
		};
		int nPoints = (int)ceil(((double)phi->GetDegree() + 1) / 2);
		return Utils::Integral(nPoints, func, x1, x2, y1, y2);
	}

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

	static double Integral(std::function<double(Point)> func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		function<double(double, double, double)> funcToIntegrate = [func](double x, double y, double z) {
			Point p(x, y, z);
			return func(p);
		};
		return Integral(funcToIntegrate, x1, x2, y1, y2, z1, z2);
	}

	static double Integral(BasisFunction<3>* phi, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		function<double(double, double, double)> func = [phi](double x, double y, double z) {
			return phi->Eval(Point(x, y, z));
		};
		int nPoints = (int)ceil(((double)phi->GetDegree() + 1) / 2);
		return Utils::Integral(nPoints, func, x1, x2, y1, y2, z1, z2);
	}

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