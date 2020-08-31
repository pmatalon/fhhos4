#pragma once
#include "../TestCase.h"
using namespace std;

class ExpSolution2DTestCase : public TestCase<2>
{
public:
	ExpSolution2DTestCase(ProblemArguments pb) :
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
		return (-pow(y, 4) - 2 * x*(1 + 2 * x*y*y))*exp(x*y*y);
	}

	static double Solution(DomPoint p)
	{
		double x = p.X;
		double y = p.Y;
		return exp(x*y*y);
	}
};
