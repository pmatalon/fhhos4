#pragma once
#include "../TestCase.h"
using namespace std;

class EDFTestCase : public TestCase<2>
{
private:
	Tensor<2>* tensorWeirdShapeInTheMiddle;
	Tensor<2>* tensorBottomStripe;
	Tensor<2>* tensorTopRectangle;
	Tensor<2>* tensorLittlePiece;
public:
	EDFTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		tensorWeirdShapeInTheMiddle = new Tensor<2>(                    1, pb.AnisotropyRatio, pb.AnisotropyAngle);
		tensorBottomStripe          = new Tensor<2>(                    7, pb.AnisotropyRatio, pb.AnisotropyAngle);
		tensorTopRectangle          = new Tensor<2>(                    4, pb.AnisotropyRatio, pb.AnisotropyAngle);
		tensorLittlePiece           = new Tensor<2>(pb.HeterogeneityRatio, pb.AnisotropyRatio, pb.AnisotropyAngle);

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
			// left
			if (pow(x-5, 2) + pow(y-10, 2) <= 8)
				return 1.0;
			// middle
			if (pow(x - 15, 2) + pow(y - 6, 2) <= 8)
				return 1.0;
			// top
			if (pow(x - 15, 2) + pow(y - 12, 2) <= 8)
				return 1.0;
			// top
			if (pow(x - 24, 2) + pow(y - 8, 2) <= 4)
				return 1.0;
			// right
			if (x >= 26 && x <= 30 && y >= 9 && y <= 13)
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

	~EDFTestCase()
	{
		delete tensorWeirdShapeInTheMiddle;
		delete tensorBottomStripe;
		delete tensorTopRectangle;
		delete tensorLittlePiece;
	}
};
