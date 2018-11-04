#pragma once
#include <functional>
#include <vector>
#include "GaussLegendre.h"

class Utils
{
public:
	// Integral on [-1, 1]
	/*static double IntegralOnReferenceInterval(std::function<double(double)> func)
	{
		return GaussLegendre(7, func);
	}*/

	// Integral on [x1, x2]
	/*static double Integral(std::function<double(double)> func, double x1, double x2)
	{
		return GaussLegendre(7, func, x1, x2);
	}*/

	// Integral on [-1, 1] x [-1, 1]
	static double IntegralOnReferenceInterval(std::function<double(double, double)> func)
	{
		GaussLegendre* gs = new GaussLegendre();
		return gs->Quadrature(func);
	}

	// Integral on [x1, x2] x [y1, y2]
	static double Integral(std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		//int nPoints = 7;
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

/*private:
	// Gauss-Legendre quadrature on [-1, 1]
	static double GaussLegendre(int nPoints, std::function<double(double)> func)
	{
		vector<double> points = GetGaussLegendrePoints(nPoints);
		vector<double> weights = GetGaussLegendreWeights(nPoints);

		double sum = 0;
		for (int i = 0; i < nPoints; i++)
			sum += func(points[i]) * weights[i];
		return sum;
	}

	// Gauss-Legendre quadrature on [x1, x2]
	static double GaussLegendre(int nPoints, std::function<double(double)> func, double x1, double x2)
	{
		vector<double> points = GetGaussLegendrePoints(nPoints);
		vector<double> weights = GetGaussLegendreWeights(nPoints);

		double sum = 0;
		for (int i = 0; i < nPoints; i++)
				sum += func((x2 - x1) / 2 * points[i] + (x1 + x2) / 2) * weights[i];
		return (x2 - x1) / 2 * sum;
	}

	// Gauss-Legendre quadrature on [-1, 1]^2
	static double GaussLegendre(int nPoints, std::function<double(double, double)> func)
	{
		vector<double> points = GetGaussLegendrePoints(nPoints);
		vector<double> weights = GetGaussLegendreWeights(nPoints);

		double sum = 0;
		for (int i = 0; i < nPoints; i++)
		{
			for (int j = 0; j < nPoints; j++)
				sum += func(points[i], points[j]) * weights[i] * weights[j];
		}
		return sum;
	}

	// Gauss-Legendre quadrature on [x1, x2] x [y1, y2]
	static double GaussLegendre(int nPoints, std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		vector<double> points = GetGaussLegendrePoints(nPoints);
		vector<double> weights = GetGaussLegendreWeights(nPoints);

		double sum = 0;
		for (int i = 0; i < nPoints; i++)
		{
			for (int j = 0; j < nPoints; j++)
				sum += func((x2 - x1) / 2 * points[i] + (x1 + x2) / 2, (y2 - y1) / 2 * points[j] + (y1 + y2) / 2) * weights[i] * weights[j];
		}
		return (x2 - x1) / 2 * (y2 - y1) / 2 * sum;
	}

	static vector<double> GetGaussLegendrePoints(int nPoints)
	{
		vector<double> points(nPoints);
		if (nPoints == 4)
		{
			points[0] = -0.3399810435848563;
			points[1] = 0.3399810435848563;
			points[2] = -0.8611363115940526;
			points[3] = 0.8611363115940526;
		}
		else if (nPoints == 5)
		{	
			points[0] = 0.0000000000000000;
			points[1] = -0.5384693101056831;
			points[2] = 0.5384693101056831;
			points[3] = -0.9061798459386640;
			points[4] = 0.9061798459386640;
		}
		else if (nPoints == 6)
		{
			points[0] = 0.6612093864662645;
			points[1] = -0.6612093864662645;
			points[2] = -0.2386191860831969;
			points[3] = 0.2386191860831969;
			points[4] = -0.9324695142031521;
			points[5] = 0.9324695142031521;
		}
		else if (nPoints == 7)
		{
			points[0] = 0.0000000000000000;
			points[1] = 0.4058451513773972;
			points[2] = -0.4058451513773972;
			points[3] = -0.7415311855993945;
			points[4] = 0.7415311855993945;
			points[5] = -0.9491079123427585;
			points[6] = 0.9491079123427585;
		}
		return points;
	}
	static vector<double> GetGaussLegendreWeights(int nPoints)
	{
		vector<double> weights(nPoints);
		if (nPoints == 4)
		{
			weights[0] = 0.6521451548625461;
			weights[1] = 0.6521451548625461;
			weights[2] = 0.3478548451374538;
			weights[3] = 0.3478548451374538;
		}
		else if (nPoints == 5)
		{
			weights[0] = 0.5688888888888889;
			weights[1] = 0.4786286704993665;
			weights[2] = 0.4786286704993665;
			weights[3] = 0.2369268850561891;
			weights[4] = 0.2369268850561891;
		}
		else if (nPoints == 6)
		{
			weights[0] = 0.3607615730481386;
			weights[1] = 0.3607615730481386;
			weights[2] = 0.4679139345726910;
			weights[3] = 0.4679139345726910;
			weights[4] = 0.1713244923791704;
			weights[5] = 0.1713244923791704;
		}
		else if (nPoints == 7)
		{	
			weights[0] = 0.4179591836734694;
			weights[1] = 0.3818300505051189;
			weights[2] = 0.3818300505051189;
			weights[3] = 0.2797053914892766;
			weights[4] = 0.2797053914892766;
			weights[5] = 0.1294849661688697;
			weights[6] = 0.1294849661688697;
		}
		return weights;
	}*/
};