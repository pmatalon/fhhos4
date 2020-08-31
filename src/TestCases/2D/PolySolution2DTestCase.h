#pragma once
#include "../TestCase.h"
using namespace std;

class PolySolution2DTestCase : public TestCase<2>
{
public:
	PolySolution2DTestCase(ProblemArguments pb) :
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

		if (pb.BCCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = BoundaryConditions::MixedConditionsExample;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}

		// Exact solution
		if (this->DiffField.IsHomogeneous && this->DiffField.IsIsotropic && pb.BCCode.compare("d") == 0)
			this->ExactSolution = this->Solution;
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
		return 2 * (y*(1 - y) + x * (1 - x));
	}

	static double Solution(DomPoint p)
	{
		double x = p.X;
		double y = p.Y;
		return x * (1 - x) * y*(1 - y);
	}
};
