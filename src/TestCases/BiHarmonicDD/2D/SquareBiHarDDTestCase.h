#pragma once
#include "../BiHarmonicDDTestCase.h"
using namespace std;

class SquareBiHarDDTestCase : public BiHarmonicDDTestCase<2>
{
public:
	SquareBiHarDDTestCase(ProblemArguments pb) :
		BiHarmonicDDTestCase()
	{
		// Boundary conditions
		this->DirichletBC = BoundaryConditions::HomogeneousDirichletEverywhere();
		this->LaplacianDirichletBC = new BoundaryConditions(BoundaryConditions::HomogeneousDirichletEverywhere());

		// Source function
		this->SourceFunction = this->SineSource2D;

		// Exact solution
		if (pb.GeoCode.compare("square") == 0)
		{
			this->MinusLaplacianOfSolution = this->SineMinusLaplacianOfSolution2D;
			this->ExactSolution = this->SineSolution2D;
		}
	}

	string Code() override
	{
		return "";
	}
	string Description() override
	{
		return "";
	}
};
