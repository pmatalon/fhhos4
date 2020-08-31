#pragma once
#include "../TestCase.h"
#include "../../Mesh/3D/CubeGeometry.h"
using namespace std;

class ExpSolution3DTestCase : public TestCase<3>
{
public:
	ExpSolution3DTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		this->DiffField = CubeGeometry::DiffField(pb.AnisotropyRatio, pb.AnisotropyAngle);

		// Source function
		this->SourceFunction = this->Source;

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0 && pb.BCCode.compare("m") != 0)
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		if (this->DiffField.IsHomogeneous && this->DiffField.IsIsotropic && pb.BCCode.compare("d") == 0)
		{
			this->ExactSolution = this->Solution;
			this->BC.DirichletFunction = this->Solution;
		}

		if (pb.BCCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = BoundaryConditions::MixedConditionsExample;
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
};
