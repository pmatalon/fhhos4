#pragma once
#include "DiffusionTestCase.h"
using namespace std;

template <int Dim>
class DefaultTestCase : public DiffusionTestCase<Dim>
{
public:
	DefaultTestCase(ProblemArguments pb) :
		DiffusionTestCase<Dim>()
	{
		// Diffusion field
		this->DiffField = DiffusionField<Dim>(pb.AnisotropyRatio, pb.AnisotropyAngle);

		// Source function
		this->SourceFunction = this->DiscontinuousSource;

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0)
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");
	}

	string Code() override
	{
		return "default";
	}
	string Description() override
	{
		return "Default";
	}
};
