#pragma once
#include "TestCase.h"
#include "../Problem/DiffusionField.h"
using namespace std;

template <int Dim>
class BiHarmonicTestCase : public TestCase<Dim>
{
public:
	DomFunction SourceFunction = nullptr;
	DomFunction MinusLaplacianOfSolution = nullptr;
	//DiffusionField<Dim> DiffField;

	BiHarmonicTestCase()
	{}

protected:
	// With homogeneous Dirichlet B.C. for both diffusion problems
	static double SineSource2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 4 * pow(4 * M_PI, 4) * sin(4 * M_PI * x) * sin(4 * M_PI * y);
	}
	static double SineMinusLaplacianOfSolution2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x) * sin(4 * M_PI * y);
	}
	static double SineSolution2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return sin(4 * M_PI * x) * sin(4 * M_PI * y);
	}

	// With mixed Dirichlet/Neumann B.C.
	static double PolySource2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 56400 * (1 - 10*x + 15*x*x) * pow(1 - y, 2) * pow(y, 4) + 18800 * x*x * (6 - 20 * x + 15 * x*x) * y*y* (6 - 20*y + 15*y*y) + 56400*pow(1-x, 2) * pow(x, 4) * (1 - 10*y + 15*y*y);
	}
	static double PolySolution2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 2350*pow(x, 4)*pow(x-1, 2)*pow(y, 4)*pow(y-1, 2);
	}

public:
	virtual ~BiHarmonicTestCase()
	{}
};
