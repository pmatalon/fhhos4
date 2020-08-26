#pragma once
#include "../TestCase.h"
using namespace std;

class SineSolution3DTestCase : public TestCase<3>
{
public:
	SineSolution3DTestCase(DiffusionPartition<3>* diffusionPartition, string bcCode) : 
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
		return "sine";
	}
	string Description() override
	{
		return "Sine function";
	}

private:
	static double Source(DomPoint p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return 3 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z);
	}

	static double Solution(DomPoint p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return sin(4 * M_PI * x)*sin(4 * M_PI * y)*sin(4 * M_PI * z);
	}

	static BoundaryConditionType MixedConditions(DomPoint p)
	{
		return p.X == 0 ? BoundaryConditionType::Dirichlet : BoundaryConditionType::Neumann;
	}
};
