#pragma once
#include "../TestCase.h"
using namespace std;

class SquareCircleTestCase : public TestCase<2>
{
private:
	Tensor<2> diffTensorSquare;
	Tensor<2> diffTensorCircle;
public:
	SquareCircleTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		diffTensorSquare = Tensor<2>(pb.HeterogeneityRatio, pb.AnisotropyRatio, pb.AnisotropyAngle);
		diffTensorCircle = Tensor<2>(1                    , pb.AnisotropyRatio, pb.AnisotropyAngle);

		this->DiffField = DiffusionField<2>("square", &diffTensorSquare, "disk", &diffTensorCircle);

		// Source function
		if (pb.SourceCode.compare("") == 0)
		{
			this->SourceFunction = [](const DomPoint& p)
			{
				double r = 6;
				double power = 1.0;
				// left
				if (Utils::IsInDisk(DomPoint(5, 10), r, p))
					return power;
				// middle
				if (Utils::IsInDisk(DomPoint(15, 6), r, p))
					return power;
				// middle
				if (Utils::IsInDisk(DomPoint(24, 8), r, p))
					return power;
				// right
				if (Utils::IsInDisk(DomPoint(28, 11), r, p))
					return power;

				// bottom
				r = 2;
				power = 7;
				if (Utils::IsInDisk(DomPoint(1, 2), r, p))
					return power;
				return 0.0;
			};
		}
		else if (pb.SourceCode.compare("sine") == 0)
			this->SourceFunction = this->SineSource2D;
		else if (pb.SourceCode.compare("poly") == 0)
			this->SourceFunction = this->PolySource2D;
		else if (pb.SourceCode.compare("exp") == 0)
			this->SourceFunction = this->ExpSource2D;
		else if (pb.SourceCode.compare("zero") == 0)
			this->SourceFunction = this->Zero;
		else
			Utils::FatalError("Unmanaged source code");

		// Boundary conditions
		if (pb.BCCode.compare("d") == 0)
		{
			// These are already the default value, but I reset them as an example of how to apply boundary conditions.
			this->BC.GetBoundaryConditionType = BoundaryConditions::DirichletEverywhere;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			/*this->BC.DirichletFunction = [&pb](const DomPoint& p)
			{
				return 1 / pb.HeterogeneityRatio;
			};*/
			this->BC.Description = "Homogeneous Dirichlet";
		}
		else if (pb.BCCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = BoundaryConditions::MixedConditionsExample;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			/*this->BC.DirichletFunction = [&pb](const DomPoint& p)
			{
				return pb.HeterogeneityRatio;
			};*/
			this->BC.NeumannFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}
		else
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		// GMSH geometric points to ignore 
		this->GeometricPointExclusionList = { 5 };
	}

	string Code() override
	{
		return "squarecircle";
	}
	string Description() override
	{
		return "Square embedding a circle";
	}
};
