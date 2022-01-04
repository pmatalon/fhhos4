#pragma once
#include "../../Utils/Utils.h"
#include "../../Mesh/PhysicalGroup.h"
using namespace std;

class BoundaryConditions
{
public:
	PbBoundaryConditions Type = PbBoundaryConditions::FullDirichlet;
	string Description;
	function<BoundaryConditionType(BoundaryGroup*)> BoundaryConditionPartition = nullptr;
	DomFunction DirichletFunction = nullptr;
	DomFunction NeumannFunction = nullptr;

	BoundaryConditions()
	{
		Type = PbBoundaryConditions::FullDirichlet;
		Description = "Dirichlet";
		BoundaryConditionPartition = DirichletEverywhere;
		DirichletFunction = Homogeneous;
		NeumannFunction = nullptr;
	}

	BoundaryConditions(PbBoundaryConditions type, function<BoundaryConditionType(BoundaryGroup*)> getBoundaryConditionType, DomFunction dirichletFunction, DomFunction neumannFunction)
	{
		this->Type = type;
		this->BoundaryConditionPartition = getBoundaryConditionType;
		this->DirichletFunction = dirichletFunction;
		this->NeumannFunction = neumannFunction;
	}

	static BoundaryConditions HomogeneousDirichletEverywhere()
	{
		BoundaryConditions bc;
		bc.Type = PbBoundaryConditions::FullDirichlet;
		bc.BoundaryConditionPartition = DirichletEverywhere;
		bc.DirichletFunction = Homogeneous;
		bc.NeumannFunction = nullptr;
		return bc;
	}

	static BoundaryConditions HomogeneousNeumannEverywhere()
	{
		BoundaryConditions bc;
		bc.Type = PbBoundaryConditions::FullNeumann;
		bc.BoundaryConditionPartition = NeumannEverywhere;
		bc.DirichletFunction = nullptr;
		bc.NeumannFunction = Homogeneous;
		return bc;
	}


	static BoundaryConditionType DirichletEverywhere(BoundaryGroup* boundaryPart)
	{
		return BoundaryConditionType::Dirichlet;
	}

	static BoundaryConditionType NeumannEverywhere(BoundaryGroup* boundaryPart)
	{
		return BoundaryConditionType::Neumann;
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