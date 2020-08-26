#pragma once
#include "../Utils/Utils.h"
#include "../Problem/DiffusionPartition.h"
#include "../Problem/BoundaryConditions.h"
using namespace std;

template <int Dim>
class TestCase
{
public:
	DomFunction SourceFunction = nullptr;
	DomFunction ExactSolution = nullptr;
	DiffusionPartition<Dim>* DiffPartition;
	BoundaryConditions BC;

	TestCase(DiffusionPartition<Dim>* diffusionPartition)
	{
		DiffPartition = diffusionPartition;
	}

	virtual string Code() = 0;
	virtual string Description() = 0;
};
