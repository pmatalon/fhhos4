#pragma once
#include "../Utils/Utils.h"
#include "../Mesh/PhysicalGroup.h"
using namespace std;

class BoundaryConditions
{
public:
	string Description;
	function<BoundaryConditionType(BoundaryGroup*)> GetBoundaryConditionType = nullptr;
	DomFunction DirichletFunction = nullptr;
	DomFunction NeumannFunction = nullptr;

	BoundaryConditions()
	{
		Description = "Dirichlet";
		GetBoundaryConditionType = DirichletEverywhere;
		DirichletFunction = Homogeneous;
		NeumannFunction = Homogeneous;
	}

	BoundaryConditions(function<BoundaryConditionType(BoundaryGroup*)> getBoundaryConditionType, DomFunction dirichletFunction, DomFunction neumannFunction)
	{
		this->GetBoundaryConditionType = getBoundaryConditionType;
		this->DirichletFunction = dirichletFunction;
		this->NeumannFunction = neumannFunction;
	}

	static BoundaryConditionType DirichletEverywhere(BoundaryGroup* boundaryPart)
	{
		return BoundaryConditionType::Dirichlet;
	}

	static BoundaryConditionType MixedConditionsExample(BoundaryGroup* boundaryPart)
	{
		return boundaryPart->Name.compare("bottomBoundary") == 0 ? BoundaryConditionType::Dirichlet : BoundaryConditionType::Neumann;
	}

	static BoundaryConditionType NeumannOnHoles(BoundaryGroup* boundaryPart)
	{
		if (boundaryPart->Name.rfind("hole", 0) == 0) // starts with "hole"
			return BoundaryConditionType::Neumann;
		return BoundaryConditionType::Dirichlet;
	}

	static double Homogeneous(const DomPoint& p)
	{
		return 0.0;
	};
};