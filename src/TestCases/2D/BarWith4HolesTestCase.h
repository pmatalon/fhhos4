#pragma once
#include "../TestCase.h"
using namespace std;

class BarWith4HolesTestCase : public TestCase<2>
{
public:
	BarWith4HolesTestCase(DiffusionField<2>* diffusionField, string bcCode) :
		TestCase(diffusionField)
	{
		this->SourceFunction = this->Source;

		if (bcCode.compare("d") == 0)
		{
			// These are already the default value, but I reset them as an example of how to apply boundary conditions.
			this->BC.GetBoundaryConditionType = BoundaryConditions::DirichletEverywhere;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Homogeneous Dirichlet";
		}
		else if (bcCode.compare("nholes") == 0)
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
