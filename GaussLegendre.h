#pragma once
#include <functional>
#include <vector>

class GaussLegendre
{
public:
	static const int MAX_POINTS = 10;
private:
	int nPoints;
	vector<double> points;
	vector<double> weights;

public:
	GaussLegendre() : GaussLegendre(MAX_POINTS)
	{ }

	GaussLegendre(int nPoints)
	{
		if (nPoints < 2)
			nPoints = 2;
		else if (nPoints > MAX_POINTS)
			nPoints = MAX_POINTS;

		this->nPoints = nPoints;

		this->points.reserve(nPoints);
		this->weights.reserve(nPoints);

		if (nPoints == 2)
		{
			points[0] = -0.5773502691896257;
			points[1] = 0.5773502691896257;

			weights[0] = 1.0000000000000000;
			weights[1] = 1.0000000000000000;
		}
		else if (nPoints == 3)
		{
			points[0] = 0.0000000000000000;
			points[1] = -0.7745966692414834;
			points[2] = 0.7745966692414834;

			weights[0] = 0.8888888888888888;
			weights[1] = 0.5555555555555556;
			weights[2] = 0.5555555555555556;
		}
		else if (nPoints == 4)
		{
			points[0] = -0.3399810435848563;
			points[1] = 0.3399810435848563;
			points[2] = -0.8611363115940526;
			points[3] = 0.8611363115940526;

			weights[0] = 0.6521451548625461;
			weights[1] = 0.6521451548625461;
			weights[2] = 0.3478548451374538;
			weights[3] = 0.3478548451374538;
		}
		else if (nPoints == 5)
		{
			points[0] = 0.0000000000000000;
			points[1] = -0.5384693101056831;
			points[2] = 0.5384693101056831;
			points[3] = -0.9061798459386640;
			points[4] = 0.9061798459386640;

			weights[0] = 0.5688888888888889;
			weights[1] = 0.4786286704993665;
			weights[2] = 0.4786286704993665;
			weights[3] = 0.2369268850561891;
			weights[4] = 0.2369268850561891;
		}
		else if (nPoints == 6)
		{
			points[0] = 0.6612093864662645;
			points[1] = -0.6612093864662645;
			points[2] = -0.2386191860831969;
			points[3] = 0.2386191860831969;
			points[4] = -0.9324695142031521;
			points[5] = 0.9324695142031521;

			weights[0] = 0.3607615730481386;
			weights[1] = 0.3607615730481386;
			weights[2] = 0.4679139345726910;
			weights[3] = 0.4679139345726910;
			weights[4] = 0.1713244923791704;
			weights[5] = 0.1713244923791704;
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

			weights[0] = 0.4179591836734694;
			weights[1] = 0.3818300505051189;
			weights[2] = 0.3818300505051189;
			weights[3] = 0.2797053914892766;
			weights[4] = 0.2797053914892766;
			weights[5] = 0.1294849661688697;
			weights[6] = 0.1294849661688697;
		}
		else if (nPoints == 8)
		{
			points[0] = -0.1834346424956498;
			points[1] = 0.1834346424956498;
			points[2] = -0.5255324099163290;
			points[3] = 0.5255324099163290;
			points[4] = -0.7966664774136267;
			points[5] = 0.7966664774136267;
			points[6] = -0.9602898564975363;
			points[7] = 0.9602898564975363;

			weights[0] = 0.3626837833783620;
			weights[1] = 0.3626837833783620;
			weights[2] = 0.3137066458778873;
			weights[3] = 0.3137066458778873;
			weights[4] = 0.2223810344533745;
			weights[5] = 0.2223810344533745;
			weights[6] = 0.1012285362903763;
			weights[7] = 0.1012285362903763;
		}
		else if (nPoints == 9)
		{
			points[0] = 0.0000000000000000;
			points[1] = -0.8360311073266358;
			points[2] = 0.8360311073266358;
			points[3] = -0.9681602395076261;
			points[4] = 0.9681602395076261;
			points[5] = -0.3242534234038089;
			points[6] = 0.3242534234038089;
			points[7] = -0.6133714327005904;
			points[8] = 0.6133714327005904;

			weights[0] = 0.3302393550012598;
			weights[1] = 0.1806481606948574;
			weights[2] = 0.1806481606948574;
			weights[3] = 0.0812743883615744;
			weights[4] = 0.0812743883615744;
			weights[5] = 0.3123470770400029;
			weights[6] = 0.3123470770400029;
			weights[7] = 0.2606106964029354;
			weights[8] = 0.2606106964029354;
		}
		else if (nPoints == 10)
		{
			points[0] = -0.1488743389816312;
			points[1] = 0.1488743389816312;
			points[2] = -0.4333953941292472;
			points[3] = 0.4333953941292472;
			points[4] = -0.6794095682990244;
			points[5] = 0.6794095682990244;
			points[6] = -0.8650633666889845;
			points[7] = 0.8650633666889845;
			points[8] = -0.9739065285171717;
			points[9] = 0.9739065285171717;

			weights[0] = 0.2955242247147529;
			weights[1] = 0.2955242247147529;
			weights[2] = 0.2692667193099963;
			weights[3] = 0.2692667193099963;
			weights[4] = 0.2190863625159820;
			weights[5] = 0.2190863625159820;
			weights[6] = 0.1494513491505806;
			weights[7] = 0.1494513491505806;
			weights[8] = 0.0666713443086881;
			weights[9] = 0.0666713443086881;
		}
	}

	/*----------*/
	/*    1D    */
	/*----------*/

	// Gauss-Legendre quadrature on [-1, 1]
	double Quadrature(std::function<double(double)> func)
	{
		double sum = 0;
		for (int i = 0; i < this->nPoints; i++)
			sum += func(this->points[i]) * this->weights[i];
		return sum;
	}

	// Gauss-Legendre quadrature on [x1, x2]
	double Quadrature(std::function<double(double)> func, double x1, double x2)
	{
		double sum = 0;
		for (int i = 0; i < this->nPoints; i++)
			sum += func((x2 - x1) / 2 * this->points[i] + (x1 + x2) / 2) * this->weights[i];
		return (x2 - x1) / 2 * sum;
	}

	/*----------*/
	/*    2D    */
	/*----------*/

	// Gauss-Legendre quadrature on [-1, 1]^2
	double Quadrature(std::function<double(double, double)> func)
	{
		double sum = 0;
		for (int i = 0; i < this->nPoints; i++)
		{
			for (int j = 0; j < this->nPoints; j++)
				sum += func(this->points[i], this->points[j]) * this->weights[i] * this->weights[j];
		}
		return sum;
	}

	// Gauss-Legendre quadrature on [x1, x2]x[y1, y2]
	double Quadrature(std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		double sum = 0;
		for (int i = 0; i < this->nPoints; i++)
		{
			for (int j = 0; j < this->nPoints; j++)
				sum += func((x2 - x1) / 2 * this->points[i] + (x1 + x2) / 2, (y2 - y1) / 2 * this->points[j] + (y1 + y2) / 2) * this->weights[i] * this->weights[j];
		}
		return (x2 - x1) / 2 * (y2 - y1) / 2 * sum;
	}
};