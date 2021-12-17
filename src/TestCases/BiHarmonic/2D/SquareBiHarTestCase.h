#pragma once
#include "../BiHarmonicTestCase.h"
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
		// Boundary conditions
		this->DirichletBC = BoundaryConditions::HomogeneousDirichletEverywhere();
		if (pb.BCCode.compare("dd") == 0)
			this->LaplacianDirichletBC = new BoundaryConditions(BoundaryConditions::HomogeneousDirichletEverywhere());
		else if (pb.BCCode.compare("dn") == 0)
			this->NeumannBC = new BoundaryConditions(BoundaryConditions::HomogeneousNeumannEverywhere());
		else
			Utils::FatalError("The requested boundary conditions (-bc " + pb.BCCode + ") are not defined for this test case.");

		// Source function
		if (pb.SourceCode.compare("") == 0)
		{
			if (pb.BCCode.compare("dd") == 0)
			{
				pb.SourceCode = "sine";
				_pb.SourceCode = "sine";
			}
			else
			{
				pb.SourceCode = "poly";
				_pb.SourceCode = "poly";
			}
		}
		if (pb.SourceCode.compare("sine") == 0)
			this->SourceFunction = this->SineSource2D;
		else if (pb.SourceCode.compare("poly") == 0)
			this->SourceFunction = this->PolySource2D;
		else
			Utils::FatalError("Unmanaged source code");

		// Exact solution
		if (pb.GeoCode.compare("square") == 0 && pb.BCCode.compare("dd") == 0 && pb.SourceCode.compare("sine") == 0)
		{
			this->MinusLaplacianOfSolution = this->SineMinusLaplacianOfSolution2D;
			this->ExactSolution = this->SineSolution2D;
		}
		else if (pb.GeoCode.compare("square") == 0 && pb.BCCode.compare("dn") == 0 && pb.SourceCode.compare("poly") == 0)
			this->ExactSolution = this->PolySolution2D;
	}

	string Code() override
	{
		if (_pb.SourceCode.compare("sine") == 0)
			return "sine";
		else if (_pb.SourceCode.compare("poly") == 0)
			return "poly";
		else
			return _pb.SourceCode;
	}
	string Description() override
	{
		if (_pb.SourceCode.compare("sine") == 0)
			return "sine solution";
		else if (_pb.SourceCode.compare("poly") == 0)
			return "polynomial solution";
		else
			return _pb.SourceCode;
	}
};
