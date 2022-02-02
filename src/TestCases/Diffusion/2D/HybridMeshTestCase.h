#pragma once
#include "../DiffusionTestCase.h"
using namespace std;

class HybridMeshTestCase : public DiffusionTestCase<2>
{
public:
	HybridMeshTestCase(ProblemArguments pb) :
		DiffusionTestCase()
	{
		// Diffusion field
		Tensor<2> tensorInterior(                    1, pb.AnisotropyRatio, pb.AnisotropyAngle); // anisotropic
		Tensor<2> tensorExterior(pb.HeterogeneityRatio,                  1,                  0); // isotropic

		this->DiffField = DiffusionField<2>("interior", tensorInterior, "exterior", tensorExterior);

		// Source function
		if (pb.SourceCode.compare("") == 0)
			pb.SourceCode = "sine";
		
		if (pb.SourceCode.compare("sine") == 0)
			this->SourceFunction = this->SineSource2D;
		else if (pb.SourceCode.compare("poly") == 0)
			this->SourceFunction = this->PolySource2D;
		else if (pb.SourceCode.compare("exp") == 0)
			this->SourceFunction = this->ExpSource2D;
		else if (pb.SourceCode.compare("one") == 0)
			this->SourceFunction = Utils::ConstantFunctionZero;
		else if (pb.SourceCode.compare("x") == 0)
			this->SourceFunction = Utils::ConstantFunctionZero;
		else if (pb.SourceCode.compare("zero") == 0)
			this->SourceFunction = Utils::ConstantFunctionZero;
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
		return "hybridmesh";
	}
	string Description() override
	{
		return "Hybrid mesh";
	}
};
