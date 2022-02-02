#pragma once
#include "../DiffusionTestCase.h"
using namespace std;

class EDFTestCase : public DiffusionTestCase<2>
{
public:
	EDFTestCase(ProblemArguments pb) :
		DiffusionTestCase()
	{
		// Diffusion field
		Tensor<2> tensorWeirdShapeInTheMiddle(                    1, pb.AnisotropyRatio, pb.AnisotropyAngle);
		Tensor<2> tensorBottomStripe         (                   30, pb.AnisotropyRatio, pb.AnisotropyAngle);
		Tensor<2> tensorTopRectangle         (pb.HeterogeneityRatio, pb.AnisotropyRatio, pb.AnisotropyAngle);
		Tensor<2> tensorLittlePiece          (                  100, pb.AnisotropyRatio, pb.AnisotropyAngle);

		map<string, Tensor<2>> tensors;
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
			this->BC = BoundaryConditions::HomogeneousDirichletEverywhere();
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
