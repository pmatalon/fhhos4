#pragma once
#include "../TestCase.h"
using namespace std;

class SquareCircleTestCase : public TestCase<2>
{
public:
	SquareCircleTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		double kappaSquare = pb.HeterogeneityRatio;
		double kappaCircle = 1;
		Tensor<2>* diffTensorSquare = new Tensor<2>(kappaSquare, pb.AnisotropyRatio, pb.AnisotropyAngle);
		Tensor<2>* diffTensorCircle = new Tensor<2>(kappaCircle, pb.AnisotropyRatio, pb.AnisotropyAngle);

		this->DiffField = DiffusionField<2>("square", diffTensorSquare, "disk", diffTensorCircle);

		// Source function
		this->SourceFunction = this->Source;

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0)
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");
	}

	string Code() override
	{
		return "squarecircle";
	}
	string Description() override
	{
		return "Square embedding a circle";
	}

private:
	static double Source(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 2 * pow(4 * M_PI, 2) * sin(4 * M_PI * x)*sin(4 * M_PI * y);
	}
};
