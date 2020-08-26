#pragma once
#include "../TestCase.h"
using namespace std;

class PolySolution1DTestCase : public TestCase<1>
{
public:
	PolySolution1DTestCase(DiffusionField<1>* diffusionField, string bcCode) :
		TestCase(diffusionField)
	{
		if (bcCode.compare("d") != 0)
			Utils::FatalError("This test case runs only with Dirichlet conditions");

		this->SourceFunction = Source;
		if (this->DiffField->IsHomogeneous)
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
