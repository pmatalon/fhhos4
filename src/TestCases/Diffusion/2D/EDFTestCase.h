#pragma once
#include "../DiffusionTestCase.h"
using namespace std;

class EDFTestCase : public DiffusionTestCase<2>
{
private:
	Tensor<2>* tensorWeirdShapeInTheMiddle;
	Tensor<2>* tensorBottomStripe;
	Tensor<2>* tensorTopRectangle;
	Tensor<2>* tensorLittlePiece;
public:
	EDFTestCase(ProblemArguments pb) :
		DiffusionTestCase()
	{
		// Diffusion field
		tensorWeirdShapeInTheMiddle = new Tensor<2>(                    1, pb.AnisotropyRatio, pb.AnisotropyAngle);
		tensorBottomStripe          = new Tensor<2>(                   30, pb.AnisotropyRatio, pb.AnisotropyAngle);
		tensorTopRectangle          = new Tensor<2>(pb.HeterogeneityRatio, pb.AnisotropyRatio, pb.AnisotropyAngle);
		tensorLittlePiece           = new Tensor<2>(                  100, pb.AnisotropyRatio, pb.AnisotropyAngle);

		map<string, Tensor<2>*> tensors;
		tensors.insert({ "bottomStripe", tensorBottomStripe });
		tensors.insert({ "topRectangle", tensorTopRectangle });
		tensors.insert({ "leftLittlePiece", tensorLittlePiece });
		tensors.insert({ "rightLittlePiece", tensorLittlePiece });
		tensors.insert({ "weirdShapeInTheMiddle", tensorWeirdShapeInTheMiddle });

		this->DiffField = DiffusionField<2>(tensors);

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
		else
			Utils::FatalError("Unmanaged source code");

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
