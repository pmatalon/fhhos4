#pragma once
#include "../TestCase.h"
using namespace std;

class HybridMeshTestCase : public TestCase<2>
{
private:
	Tensor<2>* tensorInterior;
	Tensor<2>* tensorExterior;
public:
	HybridMeshTestCase(ProblemArguments pb) :
		TestCase()
	{
		// Diffusion field
		tensorInterior = new Tensor<2>(                    1, pb.AnisotropyRatio, pb.AnisotropyAngle); // anisotropic
		tensorExterior = new Tensor<2>(pb.HeterogeneityRatio,                  1,                  0); // isotropic

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
			this->SourceFunction = this->Zero;
		else if (pb.SourceCode.compare("x") == 0)
			this->SourceFunction = this->Zero;
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
			this->BC.Description = "Homogeneous Dirichlet";
		}
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

	~HybridMeshTestCase()
	{
		delete tensorInterior;
		delete tensorExterior;
	}
};
