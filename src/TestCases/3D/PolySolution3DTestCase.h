#pragma once
#include "../TestCase.h"
#include "../../Mesh/3D/CubeGeometry.h"
using namespace std;

class PolySolution3DTestCase : public TestCase<3>
{
public:
	PolySolution3DTestCase(ProblemArguments pb) :
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
			this->ExactSolution = this->Solution;

		if (pb.BCCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = BoundaryConditions::MixedConditionsExample;
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
		double z = p.Z;
		return 2 * ((y*(1 - y)*z*(1 - z) + x * (1 - x)*z*(1 - z) + x * (1 - x)*y*(1 - y)));
	}

	static double Solution(DomPoint p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return x * (1 - x)*y*(1 - y)*z*(1 - z);
	}
};
