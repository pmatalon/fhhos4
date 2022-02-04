#pragma once
#include "../DiffusionTestCase.h"
using namespace std;

class SquareFullNeumannTestCase : public DiffusionTestCase<2>
{
private:
	ProblemArguments _pb;
public:
	SquareFullNeumannTestCase(ProblemArguments pb) :
		DiffusionTestCase(),
		_pb(pb)
	{
		// Diffusion field
		this->DiffField = DiffusionField<2>(1, 0);

		// Source function
		if (pb.SourceCode.compare("") == 0)
			pb.SourceCode = "exp";

		if (pb.SourceCode.compare("zero") == 0)
			this->SourceFunction = Utils::ConstantFunctionZero;
		else if (pb.SourceCode.compare("sine") == 0)
			this->SourceFunction = this->SineSource;
		else if (pb.SourceCode.compare("exp") == 0)
			this->SourceFunction = this->ExpSource;
		else
			Utils::FatalError("Unmanaged source code. Choices are 'zero', 'sine', 'exp'.");

		// Boundary conditions
		if (pb.BCCode.compare("d") == 0)
		{
			this->BC.Type = PbBoundaryConditions::FullDirichlet;
			this->BC.BoundaryConditionPartition = BoundaryConditions::DirichletEverywhere;
			this->BC.Description = "Dirichlet";

			if (pb.SourceCode.compare("zero") == 0)
				this->BC.DirichletFunction = this->MinusXPlusOneHalf;
			else if (pb.SourceCode.compare("sine") == 0)
				this->BC.DirichletFunction = this->SineSolution;
			else if (pb.SourceCode.compare("exp") == 0)
				this->BC.DirichletFunction = this->ExpSolution;
		}
		else if (pb.BCCode.compare("n2") == 0)
		{
			Utils::Warning("Full Neumann conditions do not yield a well-posed problem!");
			this->BC.Type = PbBoundaryConditions::FullNeumann;
			this->BC.BoundaryConditionPartition = BoundaryConditions::NeumannEverywhere;

			if (pb.SourceCode.compare("zero") == 0)
			{
				this->BC.NeumannFunction = [](const DomPoint& p)
				{
					if (abs(p.X) < Utils::NumericalZero) // x == 0
						return 1.0;
					if (abs(p.X - 1) < Utils::NumericalZero) // x == 1
						return -1.0;
					if (abs(p.Y) < Utils::NumericalZero) // y == 0
						return 0.0;
					if (abs(p.Y - 1) < Utils::NumericalZero) // x == 1
						return 0.0;
					Utils::FatalError("Not supposed to arrive.");
					return 0.0;
				};
				this->BC.Description = "Full Neumann (non homogeneous)";
			}
			else if (pb.SourceCode.compare("sine") == 0)
			{
				this->BC.NeumannFunction = Utils::ConstantFunctionZero;
				this->BC.Description = "Full Neumann (homogeneous)";
			}
			else if (pb.SourceCode.compare("exp") == 0)
			{
				this->BC.NeumannFunction = [](const DomPoint& p)
				{
					double x = p.X;
					double y = p.Y;
					if (abs(x) < Utils::NumericalZero) // x == 0
						return -y;
					if (abs(x - 1) < Utils::NumericalZero) // x == 1
						return y * exp(y);
					if (abs(y) < Utils::NumericalZero) // y == 0
						return -x;
					if (abs(y - 1) < Utils::NumericalZero) // y == 1
						return x * exp(x);
					Utils::FatalError("Not supposed to arrive.");
				};
				this->BC.Description = "Full Neumann (non-homogeneous)";
			}
		}
		else
			Utils::FatalError("The requested boundary conditions are not defined in this test case (use 'n2' or 'd').");

		// Exact solution
		if ((pb.GeoCode.compare("square") == 0 || Utils::StartsWith(pb.GeoCode, "square4quadrants")) && this->DiffField.IsHomogeneous && this->DiffField.IsIsotropic)
		{
			if (pb.SourceCode.compare("zero") == 0)
				this->ExactSolution = MinusXPlusOneHalf;
			else if (pb.SourceCode.compare("sine") == 0)
				this->ExactSolution = SineSolution;
			else if (pb.SourceCode.compare("exp") == 0)
				this->ExactSolution = ExpSolution;
		}
	}

	string Code() override
	{
		return "fullneum";
	}
	string Description() override
	{
		return "Full Neumann problem";
	}

	static double MinusXPlusOneHalf(const DomPoint& p)
	{
		return -p.X + 0.5;
	}

	static double SineSolution(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return cos(M_PI*x*x) * cos(2*M_PI * y);
	}

	static double SineSource(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 2*M_PI*cos(2*M_PI*y)*sin(M_PI*x*x) + 4*M_PI*M_PI*cos(2*M_PI*y)*cos(M_PI*x*x) + 4*x*x*M_PI*M_PI*cos(2*M_PI*y)*cos(M_PI*x*x);
	}

	static double ExpSolution(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return exp(x * y) - 1.317902151454404;
	}

	static double ExpSource(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return - x*x*exp(x*y) - y*y*exp(x*y);
	}
};
