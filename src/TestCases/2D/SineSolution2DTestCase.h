#pragma once
#include "../TestCase.h"
#include "../../Mesh/2D/SquareGeometry.h"
#include "../../Mesh/2D/Square4quadrantsGeometry.h"
using namespace std;

class SineSolution2DTestCase : public TestCase<2>
{
public:
	SineSolution2DTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		if (pb.GeoCode.compare("square") == 0)
			this->DiffField = SquareGeometry::DiffField(pb.AnisotropyRatio, pb.AnisotropyAngle);
		else if (pb.GeoCode.compare("square4quadrants") == 0)
			this->DiffField = Square4quadrantsGeometry::DiffField(pb.HeterogeneityRatio, pb.AnisotropyRatio, pb.AnisotropyAngle);
		else
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
		else if (pb.BCCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = BoundaryConditions::MixedConditionsExample;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			this->BC.NeumannFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}
		else
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		// Exact solution
		if ((pb.GeoCode.compare("square") == 0 || pb.GeoCode.compare("square4quadrants") == 0) && this->DiffField.IsHomogeneous && this->DiffField.IsIsotropic && pb.BCCode.compare("d") == 0)
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
};
