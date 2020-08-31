#pragma once
#include "../TestCase.h"
using namespace std;

class ZeroSolution2DTestCase : public TestCase<2>
{
public:
	ZeroSolution2DTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		if (pb.GeoCode.compare("square") == 0)
			this->DiffField = SquareGeometry::DiffField(pb.AnisotropyRatio, pb.AnisotropyAngle);
		else if (pb.GeoCode.compare("square4quadrants") == 0)
			this->DiffField = Square4quadrantsGeometry::DiffField(pb.HeterogeneityRatio, pb.AnisotropyRatio, pb.AnisotropyAngle);

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
		return "zero";
	}
	string Description() override
	{
		return "Zero";
	}

private:
	static double Source(DomPoint p)
	{
		return 0;
	}

	static double Solution(DomPoint p)
	{
		return 0;
	}
};
