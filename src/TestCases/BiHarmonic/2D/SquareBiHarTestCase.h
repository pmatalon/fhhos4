#pragma once
#include "../BiHarmonicTestCase.h"
using namespace std;

class SquareBiHarTestCase : public BiHarmonicTestCase<2>
{
private:
	ProblemArguments _pb;
public:
	SquareBiHarTestCase(ProblemArguments pb) :
		BiHarmonicTestCase(),
		_pb(pb)
	{
		// Boundary conditions
		this->DirichletBC = BoundaryConditions::HomogeneousDirichletEverywhere();
		this->NeumannBC = BoundaryConditions::HomogeneousNeumannEverywhere();

		// Source function
		if (pb.SourceCode.compare("") == 0)
		{
			pb.SourceCode = "sine";
			_pb.SourceCode = "sine";
		}
		if (pb.SourceCode.compare("sine") == 0)
			this->SourceFunction = this->SineSource;
		else if (pb.SourceCode.compare("poly") == 0)
			this->SourceFunction = this->PolySource;
		else
			Utils::FatalError("Unmanaged source code");

		// Exact solution
		if (pb.SourceCode.compare("sine") == 0)
			this->ExactSolution = this->SineSolution;
		else if (pb.SourceCode.compare("poly") == 0)
			this->ExactSolution = this->PolySolution;
	}

	static double SineSource(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 24 * pow(M_PI, 4) * pow(sin(M_PI * x), 2) * pow(sin(M_PI * y),2) + 8 * pow(M_PI, 4) * pow(cos(M_PI * x),2) * pow(cos(M_PI * y),2) - 16 * pow(M_PI, 4) * pow(cos(M_PI * x),2) * pow(sin(M_PI * y),2) - 16 * pow(M_PI, 4) * pow(cos(M_PI * y),2) * pow(sin(M_PI * x),2);
	}
	static double SineSolution(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return pow(sin(M_PI*x), 2) * pow(sin(M_PI*y), 2);
	}

	static double PolySource(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 56400 * (1 - 10 * x + 15 * x * x) * pow(1 - y, 2) * pow(y, 4) + 18800 * x * x * (6 - 20 * x + 15 * x * x) * y * y * (6 - 20 * y + 15 * y * y) + 56400 * pow(1 - x, 2) * pow(x, 4) * (1 - 10 * y + 15 * y * y);
		//return 18800 * pow(x,4) * pow(y,4) + 75200 * pow(x,3) * pow(y,4) * (2 * x - 2) + 112800 * pow(x,2) * pow(y,4) * pow(x - 1, 2) + 338400 * pow(x,4) * pow(y,2) * pow(x - 1, 2) + 75200 * pow(x,4) * pow(y,3) * (2 * y - 2) + 338400 * pow(x,2) * pow(y,4) * pow(y - 1, 2) + 112800 * pow(x,4) * pow(y,2) * pow(y - 1, 2) + 56400 * pow(x,4) * pow(x - 1, 2) * pow(y - 1, 2) + 56400 * pow(y,4) * pow(x - 1, 2) * pow(y - 1, 2) + 225600 * x * pow(y,4) * (2 * x - 2) * pow(y - 1, 2) + 225600 * pow(x,4) * y * (2 * y - 2) * pow(x - 1, 2) + 300800 * pow(x,3) * pow(y,3) * (2 * x - 2) * (2 * y - 2) + 451200 * pow(x,2) * pow(y,3) * (2 * y - 2) * pow(x - 1, 2) + 451200 * pow(x,3) * pow(y,2) * (2 * x - 2) * pow(y - 1, 2) + 676800 * pow(x,2) * pow(y,2) * pow(x - 1, 2) * pow(y - 1, 2);
	}
	static double PolySolution(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 2350 * pow(x, 4) * pow(x - 1, 2) * pow(y, 4) * pow(y - 1, 2);
	}

	string Code() override
	{
		if (_pb.SourceCode.compare("sine") == 0)
			return "sine";
		else if (_pb.SourceCode.compare("poly") == 0)
			return "poly";
		else
			return _pb.SourceCode;
	}
	string Description() override
	{
		if (_pb.SourceCode.compare("sine") == 0)
			return "sine solution";
		else if (_pb.SourceCode.compare("poly") == 0)
			return "polynomial solution";
		else
			return _pb.SourceCode;
	}
};
