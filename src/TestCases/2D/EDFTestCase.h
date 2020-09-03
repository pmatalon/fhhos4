#pragma once
#include "../TestCase.h"
using namespace std;

class EDFTestCase : public TestCase<2>
{
public:
	EDFTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		Tensor<2>* tensorWeirdShapeInTheMiddle = new Tensor<2>(                        1,     pb.AnisotropyRatio, pb.AnisotropyAngle);
		Tensor<2>* tensorBottomStripe          = new Tensor<2>(    pb.HeterogeneityRatio,     pb.AnisotropyRatio, pb.AnisotropyAngle);
		Tensor<2>* tensorTopRectangle          = new Tensor<2>(pow(pb.HeterogeneityRatio, 2), pb.AnisotropyRatio, pb.AnisotropyAngle);
		Tensor<2>* tensorLittlePiece           = new Tensor<2>(    pb.HeterogeneityRatio,     pb.AnisotropyRatio, pb.AnisotropyAngle);

		map<string, Tensor<2>*> tensors;
		tensors.insert({ "bottomStripe", tensorBottomStripe });
		tensors.insert({ "topRectangle", tensorTopRectangle });
		tensors.insert({ "leftLittlePiece", tensorLittlePiece });
		tensors.insert({ "rightLittlePiece", tensorLittlePiece });
		tensors.insert({ "weirdShapeInTheMiddle", tensorWeirdShapeInTheMiddle });

		this->DiffField = DiffusionField<2>(tensors);

		// Source function
		this->SourceFunction = [](DomPoint p)
		{
			double x = p.X;
			double y = p.Y;
			if (x * x + y * y <= 2)
				return 5.0;
			if (x >= 20 && x <= 25 && y >= 4 && y <= 9)
				return 1.0;
			return 0.0;
		};

		// Boundary conditions
		if (pb.BCCode.compare("d") == 0)
		{
			// These are already the default value, but I reset them as an example of how to apply boundary conditions.
			this->BC.GetBoundaryConditionType = BoundaryConditions::DirichletEverywhere;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Homogeneous Dirichlet";
		}
		else
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");
	}

	string Code() override
	{
		return "edf";
	}
	string Description() override
	{
		return "EDF";
	}
};
