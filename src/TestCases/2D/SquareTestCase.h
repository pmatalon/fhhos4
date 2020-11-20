#pragma once
#include "../TestCase.h"
#include "../../Mesher/SquareGeometry.h"
#include "../../Mesher/Square4quadrantsGeometry.h"
using namespace std;

class SquareTestCase : public TestCase<2>
{
private:
	ProblemArguments _pb;
public:
	SquareTestCase(ProblemArguments pb) :
		TestCase(),
		_pb(pb)
	{
		// Diffusion field
		if (pb.GeoCode.compare("square") == 0)
			this->DiffField = SquareGeometry::DiffField(pb.AnisotropyRatio, pb.AnisotropyAngle);
		else if (pb.GeoCode.compare("square4quadrants") == 0)
			this->DiffField = Square4quadrantsGeometry::DiffField(pb.HeterogeneityRatio, pb.AnisotropyRatio, pb.AnisotropyAngle);
		else
			this->DiffField = DiffusionField<2>(pb.AnisotropyRatio, pb.AnisotropyAngle);

		// Source function
		if (pb.SourceCode.compare("") == 0)
		{
			pb.SourceCode = "sine";
			_pb.SourceCode = "sine";
		}
		if (pb.SourceCode.compare("sine") == 0)
			this->SourceFunction = this->SineSource2D;
		else if (pb.SourceCode.compare("poly") == 0)
			this->SourceFunction = this->PolySource2D;
		else if (pb.SourceCode.compare("exp") == 0)
			this->SourceFunction = this->ExpSource2D;
		else if (pb.SourceCode.compare("one") == 0)
			this->SourceFunction = this->Zero;
		else if (pb.SourceCode.compare("x") == 0)
			this->SourceFunction = this->Zero;
		else if (pb.SourceCode.compare("zero") == 0)
			this->SourceFunction = this->Zero;
		else
			Utils::FatalError("Unmanaged source code");

		// Boundary conditions
		if (pb.BCCode.compare("d") == 0)
		{
			this->BC.GetBoundaryConditionType = BoundaryConditions::DirichletEverywhere;
			if (pb.SourceCode.compare("exp") == 0)
			{
				this->BC.DirichletFunction = this->One;
				this->BC.Description = "Dirichlet (exponential everywhere)";
			}
			else if (pb.SourceCode.compare("one") == 0)
			{
				this->BC.DirichletFunction = this->One;
				this->BC.Description = "Dirichlet (1 everywhere)";
			}
			else if (pb.SourceCode.compare("x") == 0)
			{
				this->BC.DirichletFunction = this->X;
				this->BC.Description = "Dirichlet (x everywhere)";
			}
			else
			{
				this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
				this->BC.Description = "Homogeneous Dirichlet";
			}
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
		{
			if (pb.SourceCode.compare("sine") == 0)
				this->ExactSolution = this->SineSolution2D;
			else if (pb.SourceCode.compare("poly") == 0)
				this->ExactSolution = this->PolySolution2D;
			else if (pb.SourceCode.compare("exp") == 0)
				this->ExactSolution = this->ExpSolution2D;
			else if (pb.SourceCode.compare("one") == 0)
				this->ExactSolution = this->One;
			else if (pb.SourceCode.compare("x") == 0)
				this->ExactSolution = this->X;
			else if (pb.SourceCode.compare("zero") == 0)
				this->ExactSolution = this->Zero;
		}
	}

	string Code() override
	{
		if (_pb.SourceCode.compare("sine") == 0)
			return "sine";
		else if (_pb.SourceCode.compare("poly") == 0)
			return "poly";
		else if (_pb.SourceCode.compare("exp") == 0)
			return "exp";
		else if (_pb.SourceCode.compare("one") == 0)
			return "one";
		else if (_pb.SourceCode.compare("x") == 0)
			return "x";
		else if (_pb.SourceCode.compare("zero") == 0)
			return "zero";
		else
			return _pb.SourceCode;
	}
	string Description() override
	{
		if (_pb.SourceCode.compare("sine") == 0)
			return "sine solution";
		else if (_pb.SourceCode.compare("poly") == 0)
			return "polynomial solution";
		else if (_pb.SourceCode.compare("exp") == 0)
			return "exponential solution";
		else if (_pb.SourceCode.compare("one") == 0)
			return "constant 1 solution";
		else if (_pb.SourceCode.compare("x") == 0)
			return "x solution";
		else if (_pb.SourceCode.compare("zero") == 0)
			return "zero solution";
		else
			return _pb.SourceCode;
	}
};
