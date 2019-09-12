#pragma once
#include <functional>
#include <vector>
#include <Eigen/Sparse>
#include "Types.h"
#include "GaussLegendre.h"
#include "../Mesh/Point.h"
#include "../FunctionalBasis/BasisFunction.h"

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

	static double Integral(int nPoints, std::function<double(Point)> func, double x1, double x2)
	{
		function<double(double)> funcToIntegrate = [func](double x) {
			return func(x);
		};
		return Integral(nPoints, funcToIntegrate, x1, x2);
	}

	static double Integral(std::function<double(Point)> func, double x1, double x2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2);
	}

	static double Integral(BasisFunction<1>* phi, double x1, double x2)
	{
		function<double(double)> func = [phi](double x) {
			return phi->Eval(x);
		};
		int nPoints = NumberOfRequiredQuadraturePoint(phi->GetDegree());
		return Utils::Integral(nPoints, func, x1, x2);
	}

	static double Integral(BasisFunction<1>* phi)
	{
		return Utils::Integral(phi, -1, 1);
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

	static double Integral(int nPoints, std::function<double(Point)> func, double x1, double x2, double y1, double y2)
	{
		function<double(double, double)> funcToIntegrate = [func](double x, double y) {
			Point p(x, y);
			return func(p);
		};
		return Integral(nPoints, funcToIntegrate, x1, x2, y1, y2);
	}

	static double Integral(std::function<double(Point)> func, double x1, double x2, double y1, double y2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2);
	}

	static double Integral(BasisFunction<2>* phi, double x1, double x2, double y1, double y2)
	{
		function<double(double, double)> func = [phi](double x, double y) {
			return phi->Eval(RefPoint(x, y));
		};
		int nPoints = NumberOfRequiredQuadraturePoint(phi->GetDegree());
		return Utils::Integral(nPoints, func, x1, x2, y1, y2);
	}

	static double Integral(BasisFunction<2>* phi)
	{
		return Utils::Integral(phi, -1, 1, -1, 1);
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

	static double Integral(int nPoints, std::function<double(Point)> func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		function<double(double, double, double)> funcToIntegrate = [func](double x, double y, double z) {
			Point p(x, y, z);
			return func(p);
		};
		return Integral(nPoints, funcToIntegrate, x1, x2, y1, y2, z1, z2);
	}

	static double Integral(std::function<double(Point)> func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2, z1, z2);
	}

	static double Integral(BasisFunction<3>* phi, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		function<double(double, double, double)> func = [phi](double x, double y, double z) {
			return phi->Eval(RefPoint(x, y, z));
		};
		int nPoints = NumberOfRequiredQuadraturePoint(phi->GetDegree());
		return Utils::Integral(nPoints, func, x1, x2, y1, y2, z1, z2);
	}

	static double Integral(BasisFunction<3>* phi)
	{
		return Utils::Integral(phi, -1, 1, -1, 1, -1, 1);
	}

	template <int Dim>
	static double Integral(int nPoints, std::function<double(RefPoint)> func)
	{
		if (Dim == 0)
			return func(0);
		else if (Dim == 1)
			return Integral(nPoints, func, -1, 1);
		else if (Dim == 2)
			return Integral(nPoints, func, -1, 1, -1, 1);
		else if (Dim == 3)
			return Integral(nPoints, func, -1, 1, -1, 1, -1, 1);
		else
		{
			cout << "Unmanaged dimension in Integral." << endl;
			exit(EXIT_FAILURE);
		}
	}

	template <int Dim>
	static double Integral(std::function<double(RefPoint)> func, int polynomialDegree)
	{
		int nPoints = NumberOfRequiredQuadraturePoint(polynomialDegree);
		if (Dim == 0)
			return func(0);
		else if (Dim == 1)
			return Integral(nPoints, func, -1, 1);
		else if (Dim == 2)
			return Integral(nPoints, func, -1, 1, -1, 1);
		else if (Dim == 3)
			return Integral(nPoints, func, -1, 1, -1, 1, -1, 1);
		else
		{
			cout << "Unmanaged dimension in Integral." << endl;
			exit(EXIT_FAILURE);
		}
	}

	template <int Dim>
	static double Integral(std::function<double(RefPoint)> func)
	{
		return Integral<Dim>(GaussLegendre::MAX_POINTS, func);
	}

	static int NumberOfRequiredQuadraturePoint(int polynomialDegree)
	{
		int nPoints = (int)ceil(((double)polynomialDegree + 1) / 2);
		return nPoints;
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

	static string MatrixInfo(const SparseMatrix& M, string name)
	{
		double density = (double)M.nonZeros() / (M.rows() * M.cols()) * 100;
		double roundedDensity = (int)(density * 10.0) / 10.0;
		//int roundedDensity = ceil(density);
		return "size(" + name + ")=" + to_string(M.rows()) + "x" + to_string(M.cols()) + ", \tnnz(" + name + ")=" + to_string(M.nonZeros()) + ", \tdensity(" + name + ")=" + to_string(density) + "%";
	}

};