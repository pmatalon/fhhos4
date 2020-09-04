#pragma once
#include "../Utils/Utils.h"
#include "../Problem/DiffusionField.h"
#include "../Problem/BoundaryConditions.h"
#include "../ProgramArguments.h"
using namespace std;

template <int Dim>
class TestCase
{
public:
	DomFunction SourceFunction = nullptr;
	DomFunction ExactSolution = nullptr;
	DiffusionField<Dim> DiffField;
	BoundaryConditions BC;

	TestCase()
	{}

	virtual string Code() = 0;
	virtual string Description() = 0;

	virtual ~TestCase()
	{}
};
