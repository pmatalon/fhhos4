#pragma once
#include "../DiffusionTestCase.h"
using namespace std;

class LShapeTestCase : public DiffusionTestCase<2>
{
public:
	LShapeTestCase(ProblemArguments pb) :
		DiffusionTestCase()
	{
		// Diffusion field
		if (pb.HeterogeneityRatio != 1)
			Utils::FatalError("This test case does not allow heterogeneity.");

		this->DiffField = DiffusionField<2>(pb.AnisotropyRatio, pb.AnisotropyAngle);

		// Source function
		if (pb.SourceCode.compare("") == 0 || pb.SourceCode.compare("zero") == 0)
			this->SourceFunction = Utils::ConstantFunctionZero;
		else
			Utils::FatalError("Unmanaged source code");

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0)
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		this->BC.DirichletFunction = this->Solution;
		if (this->DiffField.IsIsotropic)
			this->ExactSolution = this->Solution;
	}

	string Code() override
	{
		return "L_shape";
	}
	string Description() override
	{
		return "L-shape";
	}

	static double Solution(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;

		// Conversion to polar coordinates (r, t)
		double r = sqrt(x*x + y * y);
		double t = 0;
		if (x > 0 && y >= 0)
			t = atan(y / x);
		else if (x > 0 && y < 0)
			t = atan(y / x) + 2 * M_PI;
		else if (x < 0)
			t = atan(y / x) + M_PI;
		else if (x == 0 && y > 0)
			t = M_PI / 2;
		else if (x == 0 && y < 0)
			t = 3 * M_PI / 2;
		else
			assert(false);

		// Function in polar coordinates
		return pow(r, 2./3.) * sin(2*t/3);
	}
};
