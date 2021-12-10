#pragma once
#include "../DiffusionTestCase.h"
using namespace std;

class SquareHolesTestCase : public DiffusionTestCase<2>
{
public:
	SquareHolesTestCase(ProblemArguments pb) :
		DiffusionTestCase()
	{
		// Diffusion field
		if (pb.HeterogeneityRatio != 1)
			Utils::FatalError("This test case does not allow heterogeneity.");

		this->DiffField = DiffusionField<2>(pb.AnisotropyRatio, pb.AnisotropyAngle);

		// Source function
		if (pb.SourceCode.compare("") == 0)
			this->SourceFunction = this->PolySource2D;
		else
			Utils::FatalError("Unmanaged source code");

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
};
