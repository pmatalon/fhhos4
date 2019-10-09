#pragma once
#include "Element.h"
#include "../Problem/BoundaryConditions.h"
using namespace std;

template <int Dim>
class Face
{
private:
	BoundaryConditionType _boundaryConditionType;
	DomFunction _boundaryConditionFunction = [](DomPoint p) { return 0; };
public:
	BigNumber Number;
	bool IsDomainBoundary;
	Element<Dim>* Element1;
	Element<Dim>* Element2;

	bool IsRemovedOnCoarserGrid = false; 
	vector<Face<Dim>*> FinerFaces;

public:
	Face() { assert(false); }
	Face(BigNumber number, Element<Dim>* element1, Element<Dim>* element2)
	{
		this->Number = number;
		this->Element1 = element1;
		this->Element2 = element2;
		this->IsDomainBoundary = element2 == NULL;
		this->_boundaryConditionType = BoundaryConditionType::NotOnBoundary;
	}
	Face(BigNumber number, Element<Dim>* element1)
		:Face(number, element1, NULL)
	{
		this->IsDomainBoundary = true;
		this->_boundaryConditionType = BoundaryConditionType::Dirichlet;
	}

	void SetBoundaryConditions(BoundaryConditions* bc)
	{
		assert(IsDomainBoundary);
		if (bc->GetBoundaryConditionType(this->Center()) == BoundaryConditionType::Dirichlet)
		{
			this->_boundaryConditionType = BoundaryConditionType::Dirichlet;
			this->_boundaryConditionFunction = bc->DirichletFunction;
		}
		else
		{
			this->_boundaryConditionType = BoundaryConditionType::Neumann;
			this->_boundaryConditionFunction = bc->NeumannFunction;
		}
	}	

	bool HasDirichletBC()
	{
		return IsDomainBoundary && _boundaryConditionType == BoundaryConditionType::Dirichlet;
	}

	bool HasNeumannBC()
	{
		return IsDomainBoundary && _boundaryConditionType == BoundaryConditionType::Neumann;
	}

	bool IsBetween(Element<Dim>* element1, Element<Dim>* element2)
	{
		if (element1 == element2 && (element1 == this->Element1 || element1 == this->Element2))
			return true;
		if (element1 != this->Element1 && element1 != this->Element2)
			return false;
		if (element2 != this->Element1 && element2 != this->Element2)
			return false;
		return true;
	}

	Element<Dim>* GetNeighbour(Element<Dim>* element)
	{
		if (element == this->Element1)
			return this->Element2;
		else if (element == this->Element2)
			return this->Element1;
		return NULL;
	}

	int LocalNumberOf(Face<Dim>* finerFace)
	{
		for (int i = 0; i < this->FinerFaces.size(); i++)
		{
			if (this->FinerFaces[i] == finerFace)
				return i;
		}
		assert(false);
	}

	virtual double GetDiameter() = 0;
	virtual double Measure() = 0;
	virtual DomPoint Center() = 0;
	virtual DomPoint ConvertToDomain(RefPoint refPoint) = 0;
	virtual RefPoint ConvertToReference(DomPoint domainPoint) = 0;
	virtual double ComputeIntegral(RefFunction func) = 0;
	virtual double ComputeIntegral(RefFunction func, int polynomialDegree) = 0;
	virtual Face<Dim>* CreateSameGeometricFace(BigNumber number, Element<Dim>* element1) = 0;
	virtual void ExportFaceToMatlab(FILE* file) = 0;

	friend ostream& operator<<(ostream& os, const Face<Dim>& s)
	{
		s.Serialize(os);
		return os;
	}

	virtual void Serialize(ostream& os) const
	{
		if (Dim == 1)
			os << "Interface ";
		else if (Dim == 2)
			os << "Edge ";
		else if (Dim == 3)
			os << "Face ";
		os << this->Number << " of element";
		if (!this->IsDomainBoundary)
			os << "s";
		os << " " << this->Element1->Number;
		if (this->Element2)
			os << ", " << this->Element2->Number;
	}

	virtual ~Face() {}
};