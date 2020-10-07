#pragma once
#include "../TestCase.h"
#include "../../Mesher/SegmentGeometry.h"
using namespace std;

class Heterogeneity1DTestCase : public TestCase<1>
{
public:
	Heterogeneity1DTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		this->DiffField = SegmentGeometry::DiffField(pb.HeterogeneityRatio);

		// Source function
		this->SourceFunction = Source;

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0)
			Utils::FatalError("This test case runs only with Dirichlet conditions");

		this->ExactSolution = [this](const DomPoint& p)
		{
			double x = p.X;
			double alpha = this->DiffField.K1->LargestEigenValue;
			double a1 = -1 / (2 * alpha);
			double a2 = -0.5;
			double b1 = (1 + 3 * alpha) / (2 * alpha*(1 + alpha));
			double b2 = -(alpha + 3) / (2 * (1 + alpha));
			if (x < 0.5)
				return 4 * a1 *pow(x, 2) + 2 * b1 * x;
			else
				return 4 * a2 * pow(x - 1, 2) + 2 * b2 * (x - 1);
		};
	}

	string Code() override
	{
		return "heterog";
	}
	string Description() override
	{
		return "Heterogeneous-specific test case. The solution is a continuous piecewise polynomial function.";
	}

	static double Source(const DomPoint& p)
	{
		return 4;
	}
};
