#pragma once
#include "../DiffusionTestCase.h"
#include "../../../Mesher/SegmentGeometry.h"
using namespace std;

class SegmentTestCase : public DiffusionTestCase<1>
{
private:
	ProblemArguments _pb;
public:
	SegmentTestCase(ProblemArguments pb) :
		DiffusionTestCase(),
		_pb(pb)
	{
		// Diffusion field
		this->DiffField = SegmentGeometry::DiffField(pb.HeterogeneityRatio);

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

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0)
			Utils::FatalError("This test case runs only with Dirichlet conditions");

		// Exact solution
		if (this->DiffField.IsHomogeneous)
		{
			if (pb.SourceCode.compare("sine") == 0)
				this->ExactSolution = this->SineSolution;
			else if (pb.SourceCode.compare("poly") == 0)
				this->ExactSolution = this->PolySolution;
		}
	}

	string Code() override
	{
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

private:
	static double SineSource(const DomPoint& p)
	{
		double x = p.X;
		return sin(4 * M_PI * x);
	}
	static double SineSolution(const DomPoint& p)
	{
		double x = p.X;
		return sin(4 * M_PI * x) / (16 * pow(M_PI, 2));
	}

	static double PolySource(const DomPoint& p)
	{
		return 2;
	}
	static double PolySolution(const DomPoint& p)
	{
		double x = p.X;
		return x * (1 - x);
	}
};
