#pragma once
#include "Diffusion/BoundaryConditions.h"
#include "../ProgramArguments.h"
using namespace std;

template <int Dim>
class TestCase
{
public:
	DomFunction ExactSolution = nullptr;

	// GMSH Point ids
	vector<int> GeometricPointExclusionList;

	map<string, vector<int>> ReEntrantGeometricPoints;
	map<string, vector<string>> ReEntrantBoundary;

	TestCase()
	{}

	virtual string FilePrefix() = 0;
	virtual string Code() = 0;
	virtual string Description() = 0;
	virtual void PrintPhysicalProblem() = 0;

protected:
	static double Zero(const DomPoint& p)
	{
		return 0;
	}
	static double One(const DomPoint& p)
	{
		return 1;
	}
	static double X(const DomPoint& p)
	{
		return p.X;
	}

public:
	virtual ~TestCase()
	{}
};
