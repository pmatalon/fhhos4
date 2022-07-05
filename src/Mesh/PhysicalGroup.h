#pragma once
#include "../TestCases/Diffusion/Tensor.h"
using namespace std;

template <int Dim>
class PhysicalGroup
{
public:
	int Id;
	string Name;
	Tensor<Dim>* ConstantDiffTensor = nullptr;

	PhysicalGroup(int id)
	{
		this->Id = id;
	}

	PhysicalGroup(int id, string name) : 
		PhysicalGroup(id)
	{
		this->Name = name;
	}
};

class BoundaryGroup
{
public:
	int Id;
	string Name;
	BoundaryConditionType Condition = BoundaryConditionType::NotOnBoundary;
	DomFunction ConditionFunction = nullptr;

	BoundaryGroup(int id)
	{
		this->Id = id;
	}

	BoundaryGroup(int id, string name) :
		BoundaryGroup(id)
	{
		this->Name = name;
	}
};