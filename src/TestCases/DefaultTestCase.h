#pragma once
#include "TestCase.h"
using namespace std;

template <int Dim>
class DefaultTestCase : public TestCase<Dim>
{
public:
	DefaultTestCase(ProblemArguments pb) :
		TestCase<Dim>()
	{
		// Diffusion field
		this->DiffField = DiffusionField<Dim>(pb.AnisotropyRatio, pb.AnisotropyAngle);

		// Source function
		this->SourceFunction = [](DomPoint p)
		{
			double x = p.X;
			double y = p.Y;
			double z = p.Z;
			return x*x + y*y + z*z <= 0.1 ? 1.0 : 0.0;
		};

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
