#pragma once
#include "../TestCase.h"
#include "../../Mesh/1D/SegmentGeometry.h"
using namespace std;

class SineSolution1DTestCase : public TestCase<1>
{
public:
	SineSolution1DTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		this->DiffField = SegmentGeometry::DiffField(pb.HeterogeneityRatio);

		// Source function
		this->SourceFunction = Source;

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0)
			Utils::FatalError("This test case runs only with Dirichlet conditions");

		// Exact solution
		if (this->DiffField.IsHomogeneous)
			this->ExactSolution = Solution;
	}

	string Code() override
	{
		return "sine";
	}
	string Description() override
	{
		return "Sine function";
	}

	static double Source(DomPoint p)
	{
		double x = p.X;
		return sin(4 * M_PI * x);
	}

	static double Solution(DomPoint p)
	{
		double x = p.X;
		return sin(4 * M_PI * x) / (16 * pow(M_PI, 2));
	}
};
