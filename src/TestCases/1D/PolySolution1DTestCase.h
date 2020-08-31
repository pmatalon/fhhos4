#pragma once
#include "../TestCase.h"
#include "../../Mesh/1D/SegmentGeometry.h"
using namespace std;

class PolySolution1DTestCase : public TestCase<1>
{
public:
	PolySolution1DTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		this->DiffField = SegmentGeometry::DiffField(pb.HeterogeneityRatio);

		// Source function
		this->SourceFunction = Source;

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0)
			Utils::FatalError("This test case runs only with Dirichlet conditions");

		if (this->DiffField.IsHomogeneous)
			this->ExactSolution = Solution;
	}

	string Code() override
	{
		return "poly";
	}
	string Description() override
	{
		return "Polynomial function";
	}

	static double Source(DomPoint p)
	{
		return 2;
	}

	static double Solution(DomPoint p)
	{
		double x = p.X;
		return x * (1 - x);
	}
};
