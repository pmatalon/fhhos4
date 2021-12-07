#pragma once
#include "../../BiHarmonicTestCase.h"
#include "../../../Mesher/SquareGeometry.h"
using namespace std;

class SquareBiHarTestCase : public BiHarmonicTestCase<2>
{
private:
	ProblemArguments _pb;
public:
	SquareBiHarTestCase(ProblemArguments pb) :
		BiHarmonicTestCase(),
		_pb(pb)
	{
		// Source function
		if (pb.SourceCode.compare("") == 0)
		{
			pb.SourceCode = "sine";
			_pb.SourceCode = "sine";
		}
		if (pb.SourceCode.compare("sine") == 0)
			this->SourceFunction = this->SineSource2D;
		else
			Utils::FatalError("Unmanaged source code");

		// Boundary conditions
		if (pb.BCCode.compare("d") == 0)
		{
			this->BC.GetBoundaryConditionType = BoundaryConditions::DirichletEverywhere;
		}
		else
			Utils::FatalError("The requested boundary conditions are not defined in this test case.");

		// Exact solution
		if (pb.GeoCode.compare("square") == 0 && pb.BCCode.compare("d") == 0)
		{
			if (pb.SourceCode.compare("sine") == 0)
				this->ExactSolution = this->SineSolution2D;
		}
	}

	string Code() override
	{
		if (_pb.SourceCode.compare("sine") == 0)
			return "sine";
		else
			return _pb.SourceCode;
	}
	string Description() override
	{
		if (_pb.SourceCode.compare("sine") == 0)
			return "sine solution";
		else
			return _pb.SourceCode;
	}
};
