#pragma once
#include "../TestCase.h"
using namespace std;

class PolySolution2DTestCase : public TestCase<2>
{
public:
	PolySolution2DTestCase(DiffusionPartition<2>* diffusionPartition, string bcCode) : 
		TestCase(diffusionPartition)
	{
		if (bcCode.compare("d") != 0 && bcCode.compare("m") != 0)
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		this->SourceFunction = this->Source;
		if (this->DiffPartition->IsHomogeneous && this->DiffPartition->IsIsotropic && bcCode.compare("d") == 0)
			this->ExactSolution = this->Solution;

		if (bcCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = MixedConditions;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}
	}

	string Code() override
	{
		return "poly";
	}
	string Description() override
	{
		return "Polynomial function";
	}

private:
	static double Source(DomPoint p)
	{
		double x = p.X;
		double y = p.Y;
		return 2 * (y*(1 - y) + x * (1 - x));
	}

	static double Solution(DomPoint p)
	{
		double x = p.X;
		double y = p.Y;
		return x * (1 - x) * y*(1 - y);
	}

	static BoundaryConditionType MixedConditions(DomPoint p)
	{
		return p.X == 0 ? BoundaryConditionType::Dirichlet : BoundaryConditionType::Neumann;
	}
};
