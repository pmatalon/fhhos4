#pragma once
#include "../TestCase.h"
using namespace std;

class BarWith4HolesTestCase : public TestCase<2>
{
public:
	BarWith4HolesTestCase(ProblemArguments pb) :
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
		else if (pb.BCCode.compare("nholes") == 0)
		{
			this->BC.GetBoundaryConditionType = NeumannOnHoles;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			this->BC.NeumannFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Neumann on the holes";
		}
		else
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");
	}

	string Code() override
	{
		return "barwith4holes";
	}
	string Description() override
	{
		return "Bar with 4 holes";
	}

private:
	static double Source(DomPoint p)
	{
		double x = p.X;
		double y = p.Y;
		return 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y);
	}

	static BoundaryConditionType NeumannOnHoles(BoundaryGroup* boundaryPart)
	{
		if (boundaryPart->Name.rfind("hole", 0) == 0) // starts with "hole"
			return BoundaryConditionType::Neumann;
		return BoundaryConditionType::Dirichlet;
	}
};
