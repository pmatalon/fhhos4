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
		this->SourceFunction = [](DomPoint p)
		{
			double x = p.X;
			double y = p.Y;

			double r = 6;
			double power = 1.0;
			// left
			if (IsInDisk(DomPoint(5, 10), r, p))
				return power;
			// middle
			if (IsInDisk(DomPoint(15, 6), r, p))
				return power;
			// middle
			if (IsInDisk(DomPoint(24, 8), r, p))
				return power;
			// right
			if (IsInDisk(DomPoint(28, 11), r, p))
				return power;

			// bottom
			r = 2;
			power = 7;
			if (IsInDisk(DomPoint(1, 2), r, p))
				return power;
			/*if (IsInDisk(DomPoint(7, 2), r, p))
				return power;
			if (IsInDisk(DomPoint(13, 2), r, p))
				return power;
			if (IsInDisk(DomPoint(19, 2), r, p))
				return power;
			if (IsInDisk(DomPoint(25, 2), r, p))
				return power;
			if (IsInDisk(DomPoint(31, 2), r, p))
				return power;*/
			//if (x >= 26 && x <= 30 && y >= 9 && y <= 13)
				//return 1.0;
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

	static bool IsInDisk(DomPoint c, double r, DomPoint p)
	{
		if (pow(p.X - c.X, 2) + pow(p.Y - c.Y, 2) <= r)
			return true;
		return false;
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
