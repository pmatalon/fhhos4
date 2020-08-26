#pragma once
#include "../Utils/Utils.h"
#include "../Problem/DiffusionField.h"
#include "../Problem/BoundaryConditions.h"
using namespace std;

template <int Dim>
class TestCase
{
public:
	DomFunction SourceFunction = nullptr;
	DomFunction ExactSolution = nullptr;
	DiffusionField<Dim>* DiffField;
	BoundaryConditions BC;

	TestCase(DiffusionField<Dim>* diffusionField)
	{
		DiffField = diffusionField;
	}

	virtual string Code() = 0;
	virtual string Description() = 0;
};
