#pragma once
#include "Problem.h"
#include "../TestCases/BiHarmonicTestCase.h"
using namespace std;

template <int Dim>
class BiHarmonicProblem : public Problem<Dim>
{
protected:
	BiHarmonicTestCase<Dim>* _testCase;
	DomFunction _sourceFunction;
	BoundaryConditions* _boundaryConditions;
public:

	BiHarmonicProblem(Mesh<Dim>* mesh, BiHarmonicTestCase<Dim>* testCase, string outputDirectory)
		: Problem<Dim>(mesh, outputDirectory)
	{
		this->_testCase = testCase;
		this->_sourceFunction = testCase->SourceFunction;
		this->_boundaryConditions = &testCase->BC;

		this->_fileName = "BiHar" + to_string(Dim) + "D_" + testCase->Code();
		if (!this->_mesh->FileNamePart().empty())
			this->_fileName += "_" + this->_mesh->FileNamePart();
	}

	void PrintPhysicalProblem() override
	{
		cout << "Problem: BiHarmonic " << Dim << "D";
		cout << "    Geometry           : " << this->_mesh->GeometryDescription() << endl;
		cout << "    Test case          : " << _testCase->Description() << endl;
		cout << "    Boundary conditions: " << _testCase->BC.Description << endl;
	}

	virtual void AssertSchemeConvergence(double l2Error) {}
};