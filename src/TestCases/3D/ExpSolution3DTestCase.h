#pragma once
#include "../TestCase.h"
using namespace std;

class ExpSolution3DTestCase : public TestCase<3>
{
public:
	ExpSolution3DTestCase(DiffusionField<3>* diffusionField, string bcCode) : 
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
			this->BC.GetBoundaryConditionType = MixedConditions;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}
	}

	string Code() override
	{
		return "exp";
	}
	string Description() override
	{
		return "Exponential function";
	}

private:
	static double Source(DomPoint p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return -(pow(y, 4)*pow(z, 6) + 2 * x*pow(z, 3) + 4 * x*x*y*y*pow(z, 6) + 6 * x*y*y*z + 9 * x*x*pow(y, 4)*pow(z, 4))*exp(x*y*y*z*z*z);
	}

	static double Solution(DomPoint p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return exp(x*y*y*z*z*z);
	}

	static BoundaryConditionType MixedConditions(DomPoint p)
	{
		return p.X == 0 ? BoundaryConditionType::Dirichlet : BoundaryConditionType::Neumann;
	}
};
