#pragma once
#include <string>
#include <vector>
#include "../Utils/Types.h"
using namespace std;

class PhysicalGroup
{
public:
	int Id;
	string Name;

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

class BoundaryGroup : public PhysicalGroup
{
public:
	BoundaryGroup(int id) : PhysicalGroup(id) {}
	BoundaryGroup(int id, string name) : PhysicalGroup(id, name) {}

	BoundaryConditionType Condition = BoundaryConditionType::NotOnBoundary;
	DomFunction ConditionFunction = nullptr;
};