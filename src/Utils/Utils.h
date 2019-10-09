#pragma once
#include <functional>
#include <vector>
#include "Types.h"
#include "GaussLegendre.h"
#include "../Mesh/Point.h"
#include "../FunctionalBasis/BasisFunction.h"
using namespace std;

using RefFunction = function<double(RefPoint)>;
using DomFunction = function<double(DomPoint)>;

class Utils
{
public:
	//-------------//
	// Integral 1D //
	//-------------//

	static double Integral(int nPoints, std::function<double(double)> func, double x1, double x2)
	{
		GaussLegendre gs(nPoints);
		return gs.Quadrature(func, x1, x2);
	}
	
	static double Integral(std::function<double(double)> func, double x1, double x2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2);
	}

	static double Integral(int nPoints, DomFunction func, double x1, double x2)
	{
		function<double(double)> funcToIntegrate = [func](double x) {
			return func(DomPoint(x));
		};
		return Integral(nPoints, funcToIntegrate, x1, x2);
	}

	static double Integral(DomFunction func, double x1, double x2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2);
	}

	//-------------//
	// Integral 2D //
	//-------------//

	// Integral on [x1, x2] x [y1, y2]
	static double Integral(int nPoints, std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		GaussLegendre gs(nPoints);
		return gs.Quadrature(func, x1, x2, y1, y2);
	}

	static double Integral(std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2);
	}

	static double Integral(int nPoints, DomFunction func, double x1, double x2, double y1, double y2)
	{
		function<double(double, double)> funcToIntegrate = [func](double x, double y) {
			return func(DomPoint(x, y));
		};
		return Integral(nPoints, funcToIntegrate, x1, x2, y1, y2);
	}

	static double Integral(DomFunction func, double x1, double x2, double y1, double y2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2);
	}

	//-------------//
	// Integral 3D //
	//-------------//

	// Integral on [x1, x2] x [y1, y2]
	static double Integral(int nPoints, std::function<double(double, double, double)> func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		GaussLegendre gs(nPoints);
		return gs.Quadrature(func, x1, x2, y1, y2, z1, z2);
	}

	static double Integral(std::function<double(double, double, double)> func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2, z1, z2);
	}

	static double Integral(int nPoints, DomFunction func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		function<double(double, double, double)> funcToIntegrate = [func](double x, double y, double z) {
			return func(DomPoint(x, y, z));
		};
		return Integral(nPoints, funcToIntegrate, x1, x2, y1, y2, z1, z2);
	}

	static double Integral(DomFunction func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2, z1, z2);
	}


	template <int Dim>
	static double Integral(int nPoints, RefFunction func)
	{
		if (Dim == 0)
			return func(0);
		GaussLegendre gs(nPoints);
		return gs.QuadratureDim<Dim>(func);
	}

	template <int Dim>
	static double Integral(RefFunction func)
	{
		return Integral<Dim>(GaussLegendre::MAX_POINTS, func);
	}

	template <int Dim>
	static double Integral(BasisFunction<Dim>* phi)
	{
		RefFunction func = [phi](RefPoint p) {
			return phi->Eval(p);
		};
		int nPoints = NumberOfRequiredQuadraturePoint(phi->GetDegree());
		return Utils::Integral<Dim>(nPoints, func);
	}

	template <int Dim>
	static double Integral(RefFunction func, int polynomialDegree)
	{
		if (Dim == 0)
			return func(0);
		int nPoints = NumberOfRequiredQuadraturePoint(polynomialDegree);
		return Integral<Dim>(nPoints, func);
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