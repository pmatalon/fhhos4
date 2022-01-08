#pragma once
#include "../DiffusionTestCase.h"
#include "../../../Mesher/CubeGeometry.h"
using namespace std;

class CubeTestCase : public DiffusionTestCase<3>
{
private:
	ProblemArguments _pb;
public:
	CubeTestCase(ProblemArguments pb) :
		DiffusionTestCase(),
		_pb(pb)
	{
		// Diffusion field
		this->DiffField = CubeGeometry::DiffField(pb.AnisotropyRatio, pb.AnisotropyAngle);

		// Source function
		if (pb.SourceCode.compare("") == 0)
		{
			pb.SourceCode = "sine";
			_pb.SourceCode = "sine";
		}
		if (pb.SourceCode.compare("sine") == 0)
			this->SourceFunction = this->SineSource3D;
		else if (pb.SourceCode.compare("poly") == 0)
			this->SourceFunction = this->PolySource3D;
		else if (pb.SourceCode.compare("exp") == 0)
			this->SourceFunction = this->ExpSource3D;
		else if (pb.SourceCode.compare("zero") == 0)
			this->SourceFunction = Utils::ConstantFunctionZero;
		else
			Utils::FatalError("Unmanaged source code");

		// Boundary conditions
		if (pb.BCCode.compare("d") != 0 && pb.BCCode.compare("m") != 0)
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		if (pb.BCCode.compare("m") == 0)
		{
			this->BC.Type = PbBoundaryConditions::MixedDirichletNeumann;
			this->BC.BoundaryConditionPartition = BoundaryConditions::MixedConditionsExample;
			this->BC.Description = "Mixed Neumann-Dirichlet";
		}

		if (pb.GeoCode.compare("cube") == 0 && this->DiffField.IsHomogeneous && this->DiffField.IsIsotropic && pb.BCCode.compare("d") == 0)
		{
			if (pb.SourceCode.compare("sine") == 0)
				this->ExactSolution = this->SineSolution3D;
			else if (pb.SourceCode.compare("poly") == 0)
				this->ExactSolution = this->PolySolution3D;
			else if (pb.SourceCode.compare("exp") == 0)
				this->ExactSolution = this->ExpSolution3D;
			else if (pb.SourceCode.compare("zero") == 0)
				this->ExactSolution = Utils::ConstantFunctionZero;
		}
	}

	string Code() override
	{
		return _pb.SourceCode;
	}
	string Description() override
	{
		if (_pb.SourceCode.compare("sine") == 0)
			return "sine solution";
		else if (_pb.SourceCode.compare("poly") == 0)
			return "polynomial solution";
		else if (_pb.SourceCode.compare("exp") == 0)
			return "exponential solution";
		else if (_pb.SourceCode.compare("zero") == 0)
			return "zero solution";
		else
			return _pb.SourceCode;
	}
};
