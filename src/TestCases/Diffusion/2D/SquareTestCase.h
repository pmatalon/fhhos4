#pragma once
#include "../DiffusionTestCase.h"
#include "../../../Mesher/SquareGeometry.h"
#include "../../../Mesher/Square4quadrantsGeometry.h"
using namespace std;

class SquareTestCase : public DiffusionTestCase<2>
{
private:
	ProblemArguments _pb;
public:
	SquareTestCase(ProblemArguments pb) :
		DiffusionTestCase(),
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
			this->SourceFunction = Utils::ConstantFunctionZero;
		else if (pb.SourceCode.compare("x") == 0)
			this->SourceFunction = Utils::ConstantFunctionZero;
		else if (pb.SourceCode.compare("zero") == 0)
			this->SourceFunction = Utils::ConstantFunctionZero;
		else if (pb.SourceCode.compare("radiator") == 0)
			this->SourceFunction = RadiatorSource;
		else if (pb.SourceCode.compare("tworadiators") == 0)
			this->SourceFunction = TwoRadiatorSource;
		else
			Utils::FatalError("Unmanaged source code");

		// Boundary conditions
		if (pb.BCCode.compare("d") == 0)
		{
			this->BC.Type = PbBoundaryConditions::FullDirichlet;
			this->BC.BoundaryConditionPartition = BoundaryConditions::DirichletEverywhere;
			if (pb.SourceCode.compare("exp") == 0)
			{
				this->BC.DirichletFunction = this->ExpSolution2D;
				this->BC.Description = "Dirichlet (exponential everywhere)";
			}
			else if (pb.SourceCode.compare("one") == 0)
			{
				this->BC.DirichletFunction = Utils::ConstantFunctionOne;
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
			this->BC.Type = PbBoundaryConditions::MixedDirichletNeumann;
			this->BC.BoundaryConditionPartition = BoundaryConditions::MixedConditionsExample;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			this->BC.NeumannFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}
		else if (pb.BCCode.compare("n") == 0)
		{
			Utils::Warning("Full Neumann conditions do not yield a well-posed problem!");
			this->BC.Type = PbBoundaryConditions::FullNeumann;
			this->BC.BoundaryConditionPartition = BoundaryConditions::NeumannEverywhere;
			this->BC.NeumannFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Full Neumann (homogeneous)";
		}
		else
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		// Exact solution
		if ((pb.GeoCode.compare("square") == 0 || Utils::StartsWith(pb.GeoCode, "square4quadrants")) && this->DiffField.IsHomogeneous && this->DiffField.IsIsotropic)
		{
			if (pb.BCCode.compare("d") == 0)
			{
				if (pb.SourceCode.compare("sine") == 0)
				{
					this->ExactSolution = this->SineSolution2D;
					this->ExactSolution_Neumann = this->SineSolution2D_Neumann;
				}
				else if (pb.SourceCode.compare("poly") == 0)
					this->ExactSolution = this->PolySolution2D;
				else if (pb.SourceCode.compare("exp") == 0)
					this->ExactSolution = this->ExpSolution2D;
				else if (pb.SourceCode.compare("one") == 0)
					this->ExactSolution = Utils::ConstantFunctionOne;
				else if (pb.SourceCode.compare("x") == 0)
					this->ExactSolution = this->X;
				else if (pb.SourceCode.compare("zero") == 0)
					this->ExactSolution = Utils::ConstantFunctionZero;
			}
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

	static double RadiatorSource(const DomPoint& p)
	{
		if (p.X >= 0.6 && p.X <= 0.9 && p.Y >= 0.2 && p.Y <= 0.3)
			return 20;
		return 0;
	}

	static double TwoRadiatorSource(const DomPoint& p)
	{
		if (p.X >= 0.6 && p.X <= 0.9 && p.Y >= 0.2 && p.Y <= 0.3)
			return 10;
		if (p.X >= 0.05 && p.X <= 0.15 && p.Y >= 0.4 && p.Y <= 0.7)
			return 10;
		return 0;
	}
};
