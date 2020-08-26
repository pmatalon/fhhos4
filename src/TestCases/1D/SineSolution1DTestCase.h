#pragma once
#include "../TestCase.h"
using namespace std;

class SineSolution1DTestCase : public TestCase<1>
{
public:
	SineSolution1DTestCase(DiffusionPartition<1>* diffusionPartition, string bcCode) :
		TestCase(diffusionPartition)
	{
		if (bcCode.compare("d") != 0)
			Utils::FatalError("This test case runs only with Dirichlet conditions");

		this->SourceFunction = Source;
		if (this->DiffPartition->IsHomogeneous)
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
