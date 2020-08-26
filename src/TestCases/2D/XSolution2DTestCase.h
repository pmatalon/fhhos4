#pragma once
#include "../TestCase.h"
using namespace std;

class XSolution2DTestCase : public TestCase<2>
{
public:
	XSolution2DTestCase(DiffusionPartition<2>* diffusionPartition, string bcCode) : 
		TestCase(diffusionPartition)
	{
		if (bcCode.compare("d") != 0 && bcCode.compare("m") != 0)
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		this->SourceFunction = this->Source;
		if (this->DiffPartition->IsHomogeneous && this->DiffPartition->IsIsotropic && bcCode.compare("d") == 0)
		{
			this->ExactSolution = this->Solution;
			this->BC.DirichletFunction = this->Solution;
		}

		if (bcCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = MixedConditions;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}
	}

	string Code() override
	{
		return "x";
	}
	string Description() override
	{
		return "x";
	}

private:
	static double Source(DomPoint p)
	{
		return 0;
	}

	static double Solution(DomPoint p)
	{
		return p.X;
	}

	static BoundaryConditionType MixedConditions(DomPoint p)
	{
		return p.X == 0 ? BoundaryConditionType::Dirichlet : BoundaryConditionType::Neumann;
	}
};
