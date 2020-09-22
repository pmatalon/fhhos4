#pragma once
#include "../TestCase.h"
using namespace std;

class SquareHolesTestCase : public TestCase<2>
{
public:
	SquareHolesTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		if (pb.HeterogeneityRatio != 1)
			Utils::FatalError("This test case does not allow heterogeneity.");

		this->DiffField = DiffusionField<2>(pb.AnisotropyRatio, pb.AnisotropyAngle);

		// Source function
		this->SourceFunction = this->Source;

		// Boundary conditions
		if (pb.BCCode.compare("d") == 0)
		{
			// These are already the default value, but I reset them as an example of how to apply boundary conditions.
			this->BC.GetBoundaryConditionType = BoundaryConditions::DirichletEverywhere;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Homogeneous Dirichlet";
		}
		else if (pb.BCCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = BoundaryConditions::NeumannOnHoles;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			this->BC.NeumannFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Neumann on the holes";
		}
		else
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");
	}

	string Code() override
	{
		return "squareholes";
	}
	string Description() override
	{
		return "Rectangle with 24 square holes";
	}

private:
	static double Source(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 2 * (y*(1 - y) + x * (1 - x));
	}
};
