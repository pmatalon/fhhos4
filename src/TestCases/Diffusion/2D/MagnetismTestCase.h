#pragma once
#include "../DiffusionTestCase.h"
using namespace std;

class MagnetismTestCase : public DiffusionTestCase<2>
{
private:
	Tensor<2> tensorExterior;
	Tensor<2> tensorMiddle;
	Tensor<2> tensorInterior;
	Tensor<2> tensorLittlePieces;
public:
	MagnetismTestCase(ProblemArguments pb) :
		DiffusionTestCase()
	{
		// Diffusion field
		tensorExterior     = Tensor<2>(pb.HeterogeneityRatio      , pb.AnisotropyRatio, pb.AnisotropyAngle);
		tensorMiddle       = Tensor<2>(                          1, pb.AnisotropyRatio, pb.AnisotropyAngle);
		tensorInterior     = Tensor<2>(pb.HeterogeneityRatio      , pb.AnisotropyRatio, pb.AnisotropyAngle);
		tensorLittlePieces = Tensor<2>(sqrt(pb.HeterogeneityRatio), pb.AnisotropyRatio, pb.AnisotropyAngle);

		map<string, Tensor<2>*> tensors;
		tensors.insert({ "Exterior", &tensorExterior });
		tensors.insert({ "Middle", &tensorMiddle });
		tensors.insert({ "Interior", &tensorInterior });
		tensors.insert({ "LittlePieces", &tensorLittlePieces });

		this->DiffField = DiffusionField<2>(tensors);

		// Source function
		if (pb.SourceCode.compare("") == 0)
		{
			this->SourceFunction = [](const DomPoint& p)
			{
				double r = 6;
				double power = 1.0;
				/*// left
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
					return power;*/

					// bottom
				/*r = 2;
				power = 7;
				if (Utils::IsInDisk(DomPoint(1, 2), r, p))*/
				r = 0.01;
				power = 1;// 1e-6;
				if (Utils::IsInDisk(DomPoint(0, 0.8), r, p))
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
			this->BC.GetBoundaryConditionType = BoundaryConditions::DirichletEverywhere;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;//this->SineSolution2D; //BoundaryConditions::Homogeneous;
			this->BC.Description = "Homogeneous Dirichlet";
		}
		else if (pb.BCCode.compare("m") == 0)
		{
			this->BC.GetBoundaryConditionType = BoundaryConditions::MixedConditionsExample;
			this->BC.DirichletFunction = BoundaryConditions::Homogeneous;
			this->BC.NeumannFunction = BoundaryConditions::Homogeneous;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}
		else
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		// GMSH geometric points to ignore 
		this->GeometricPointExclusionList = { 1, 2, 3, 6, 37, 39, 100 };

		// Re-entrant corners
		this->ReEntrantGeometricPoints.insert({ "Middle", {24, 34} });
		this->ReEntrantGeometricPoints.insert({ "Exterior", {21, 31} });
	}

	string Code() override
	{
		return "magnetism";
	}
	string Description() override
	{
		return "Magnetism";
	}
};
