#pragma once
#include "BiHarmonicTestCase.h"
using namespace std;

template <int Dim>
class BiHarDefaultTestCase : public BiHarmonicTestCase<Dim>
{
public:
	BiHarDefaultTestCase(ProblemArguments pb)
	{
		// Boundary conditions
		this->DirichletBC = BoundaryConditions::HomogeneousDirichletEverywhere();
		this->NeumannBC = BoundaryConditions::HomogeneousNeumannEverywhere();

		// Source function
		this->SourceFunction = Utils::ConstantFunctionOne;

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0)
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");
	}

	string Code() override
	{
		return "default";
	}
	string Description() override
	{
		return "Default";
	}
};
