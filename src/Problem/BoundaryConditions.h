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
	string Description;
	function<BoundaryConditionType(DomPoint)> GetBoundaryConditionType = nullptr;
	DomFunction DirichletFunction = nullptr;
	DomFunction NeumannFunction = nullptr;

	BoundaryConditions()
	{
		Description = "Dirichlet";
		GetBoundaryConditionType = DirichletEverywhere;
		DirichletFunction = Homogeneous;
		NeumannFunction = Homogeneous;
	}

	BoundaryConditions(function<BoundaryConditionType(DomPoint)> getBoundaryConditionType, DomFunction dirichletFunction, DomFunction neumannFunction)
	{
		this->GetBoundaryConditionType = getBoundaryConditionType;
		this->DirichletFunction = dirichletFunction;
		this->NeumannFunction = neumannFunction;
	}

	static BoundaryConditionType DirichletEverywhere(DomPoint p)
	{
		return BoundaryConditionType::Dirichlet;
	};

	static double Homogeneous(DomPoint p)
	{
		return 0.0;
	};
};