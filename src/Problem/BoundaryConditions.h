#pragma once
#include <functional>
#include "../Mesh/Point.h"
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
	function<double(DomPoint)> DirichletFunction;
	function<double(DomPoint)> NeumannFunction;
	bool HomogeneousDirichlet = false;

	BoundaryConditions(function<BoundaryConditionType(DomPoint)> getBoundaryConditionType, function<double(DomPoint)> dirichletFunction, function<double(DomPoint)> neumannFunction)
	{
		this->GetBoundaryConditionType = getBoundaryConditionType;
		this->DirichletFunction = dirichletFunction;
		this->NeumannFunction = neumannFunction;
	}
};