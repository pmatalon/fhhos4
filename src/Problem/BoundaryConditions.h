#pragma once
#include "../Utils/Utils.h"
using namespace std;

enum class BoundaryConditionType : unsigned
{
	NotOnBoundary = 0,
	Dirichlet = 1,
	Neumann = 2
};

class BoundaryConditions
{
public:
	function<BoundaryConditionType(DomPoint)> GetBoundaryConditionType;
	DomFunction DirichletFunction;
	DomFunction NeumannFunction;
	bool HomogeneousDirichlet = false;

	BoundaryConditions(function<BoundaryConditionType(DomPoint)> getBoundaryConditionType, DomFunction dirichletFunction, DomFunction neumannFunction)
	{
		this->GetBoundaryConditionType = getBoundaryConditionType;
		this->DirichletFunction = dirichletFunction;
		this->NeumannFunction = neumannFunction;
	}
};