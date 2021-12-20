#pragma once
#include "../BiHarmonicTestCase.h"
using namespace std;

class SquareBiHarTestCase : public BiHarmonicTestCase<2>
{
public:
	SquareBiHarTestCase(ProblemArguments pb) :
		BiHarmonicTestCase()
	{
		// Boundary conditions
		this->DirichletBC = BoundaryConditions::HomogeneousDirichletEverywhere();
		this->NeumannBC = BoundaryConditions::HomogeneousNeumannEverywhere();

		// Source function
		this->SourceFunction = this->PolySource2D;

		// Exact solution
		this->ExactSolution = this->PolySolution2D;
	}

	// With mixed Dirichlet/Neumann B.C.
	static double PolySource2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 56400 * (1 - 10 * x + 15 * x * x) * pow(1 - y, 2) * pow(y, 4) + 18800 * x * x * (6 - 20 * x + 15 * x * x) * y * y * (6 - 20 * y + 15 * y * y) + 56400 * pow(1 - x, 2) * pow(x, 4) * (1 - 10 * y + 15 * y * y);
	}
	static double PolySolution2D(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 2350 * pow(x, 4) * pow(x - 1, 2) * pow(y, 4) * pow(y - 1, 2);
	}

	string Code() override
	{
		return "poly";
	}
	string Description() override
	{
		return "polynomial solution";
	}
};
