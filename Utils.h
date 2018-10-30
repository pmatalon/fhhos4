#pragma once
#include <functional>


class Utils
{
public:
	static double Integral(std::function<double(double)> func, double x1, double x2)
	{
		return GaussLegendre4Point(func, x1, x2);
	}

	static double Integral(std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		if (x1 == x2)
		{
			std::function<double(double)> func1D = [func, x1](double y) {
				return func(x1, y);
			};
			return GaussLegendre4Point(func1D, y1, y2);
		}
		if (y1 == y2)
		{
			std::function<double(double)> func1D = [func, y1](double x) {
				return func(x, y1);
			};
			return GaussLegendre4Point(func1D, x1, x2);
		}
		return GaussLegendre4Point(func, x1, x2, y1, y2);
	}

private:
	static double GaussLegendre4Point(std::function<double(double)> func, double x1, double x2)
	{
		// Gauss-Legendre quadrature (4 points) on [-1, 1]

		double points[4];
		double weights[4];

		points[0] = -0.3399810435848563;
		points[1] =  0.3399810435848563;
		points[2] = -0.8611363115940526;
		points[3] =  0.8611363115940526;

		weights[0] = 0.6521451548625461;
		weights[1] = 0.6521451548625461;
		weights[2] = 0.3478548451374538;
		weights[3] = 0.3478548451374538;

		double sum = 0;
		for (int i = 0; i < 4; i++)
				sum += func((x2 - x1) / 2 * points[i] + (x1 + x2) / 2) * weights[i];
		return (x2 - x1) / 2 * sum;
	}

	static double GaussLegendre4Point(std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		// Gauss-Legendre quadrature (4 points) on [-1, 1]^2

		double points[4];
		double weights[4];

		points[0] = -0.3399810435848563;
		points[1] = 0.3399810435848563;
		points[2] = -0.8611363115940526;
		points[3] = 0.8611363115940526;

		weights[0] = 0.6521451548625461;
		weights[1] = 0.6521451548625461;
		weights[2] = 0.3478548451374538;
		weights[3] = 0.3478548451374538;

		double sum = 0;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
				sum += func((x2 - x1) / 2 * points[i] + (x1 + x2) / 2, (y2 - y1) / 2 * points[j] + (y1 + y2) / 2) * weights[i] * weights[j];
		}
		return (x2 - x1) / 2 * (y2 - y1) / 2 * sum;
	}
};