#pragma once
#include "../TestCase.h"
#include "../../Mesh/3D/CubeGeometry.h"
using namespace std;

class SineSolution3DTestCase : public TestCase<3>
{
public:
	SineSolution3DTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		this->DiffField = CubeGeometry::DiffField(pb.AnisotropyRatio, pb.AnisotropyAngle);

		// Source function
		this->SourceFunction = this->Source;

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0 && pb.BCCode.compare("m") != 0)
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		if (pb.BCCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = BoundaryConditions::MixedConditionsExample;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}

		if (pb.GeoCode.compare("cube") == 0 && this->DiffField.IsHomogeneous && this->DiffField.IsIsotropic && pb.BCCode.compare("d") == 0)
			this->ExactSolution = this->Solution;
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
};
