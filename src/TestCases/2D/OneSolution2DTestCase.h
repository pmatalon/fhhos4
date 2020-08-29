#pragma once
#include "../TestCase.h"
using namespace std;

class OneSolution2DTestCase : public TestCase<2>
{
public:
	OneSolution2DTestCase(DiffusionField<2>* diffusionField, string bcCode) : 
		TestCase(diffusionField)
	{
		if (bcCode.compare("d") != 0 && bcCode.compare("m") != 0)
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		this->SourceFunction = this->Source;
		if (this->DiffField->IsHomogeneous && this->DiffField->IsIsotropic && bcCode.compare("d") == 0)
		{
			this->ExactSolution = this->Solution;
			this->BC.DirichletFunction = this->Solution;
		}

		if (bcCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = BoundaryConditions::MixedConditionsExample;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}
	}

	string Code() override
	{
		return "one";
	}
	string Description() override
	{
		return "Constant 1";
	}

private:
	static double Source(DomPoint p)
	{
		return 0;
	}

	static double Solution(DomPoint p)
	{
		return 1;
	}
};
