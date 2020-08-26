#pragma once
#include "../TestCase.h"
using namespace std;

class SineSolution2DTestCase : public TestCase<2>
{
public:
	SineSolution2DTestCase(DiffusionField<2>* diffusionField, string bcCode) : 
		TestCase(diffusionField)
	{
		this->SourceFunction = this->Source;

		if (this->DiffField->IsHomogeneous && this->DiffField->IsIsotropic && bcCode.compare("d") == 0)
			this->ExactSolution = this->Solution;

		if (bcCode.compare("d") == 0)
		{
			// These are already the default value, but I reset them as an example of how to apply boundary conditions.
			this->BC.GetBoundaryConditionType = BoundaryConditions::DirichletEverywhere;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Homogeneous Dirichlet";
		}
		else if (bcCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = MixedConditions;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			this->BC.NeumannFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}
		else
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");
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
		return 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y);
	}

	static double Solution(DomPoint p)
	{
		double x = p.X;
		double y = p.Y;
		double a = 1;//anisotropyCoefficients1[0];
		double b = 1;//anisotropyCoefficients1[1];
		return 2 / (a + b) * sin(4 * M_PI * x)*sin(4 * M_PI * y);
	}

	static BoundaryConditionType MixedConditions(DomPoint p)
	{
		return p.X == 0 ? BoundaryConditionType::Dirichlet : BoundaryConditionType::Neumann;
	}
};
