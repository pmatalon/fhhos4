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
	Face<Dim>* CoarseFace = nullptr;

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

	//-----------------------//
	//   Virtual functions   //
	//-----------------------//

	virtual GeometricShapeWithReferenceShape<Dim-1>* Shape() const = 0;

	virtual void ExportFaceToMatlab(FILE* file) = 0;

	//---------------------------//
	//   Geometric information   //
	//---------------------------//

	// Geometric information
	virtual double Diameter() const
	{
		return Shape()->Diameter();
	}
	virtual double Measure() const
	{
		return Shape()->Measure();
	}
	virtual DomPoint Center() const
	{
		return Shape()->Center();
	}
	virtual bool Contains(DomPoint p) const
	{
		return Shape()->Contains(p);
	}

	bool HasVertex(Vertex* v, bool compareCoordinates = false)
	{
		return this->Shape()->HasVertex(v, compareCoordinates);
	}

	bool HasSameVertices(Face<Dim>* other, bool compareCoordinates = false)
	{
		return this->Shape()->HasSameVertices(other->Shape(), compareCoordinates);
	}

	virtual Face<Dim>* CreateSameGeometricFace(BigNumber number, Element<Dim>* element1) = 0;

	// Transformation to reference element
	virtual DomPoint ConvertToDomain(RefPoint refPoint) const
	{
		return Shape()->ConvertToDomain(refPoint);
	}
	virtual RefPoint ConvertToReference(DomPoint domainPoint) const
	{
		return Shape()->ConvertToReference(domainPoint);
	}

	// Integral
	virtual double Integral(RefFunction func) const
	{
		return Shape()->Integral(func);
	}
	virtual double Integral(RefFunction func, int polynomialDegree) const
	{
		return Shape()->Integral(func, polynomialDegree);
	}

	//---------------------------------------//
	//   Relations to other faces/elements   //
	//---------------------------------------//

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

	//-----------------------------//
	//     Problem information     //
	//-----------------------------//

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

	inline bool HasDirichletBC()
	{
		return IsDomainBoundary && _boundaryConditionType == BoundaryConditionType::Dirichlet;
	}

	inline bool HasNeumannBC()
	{
		return IsDomainBoundary && _boundaryConditionType == BoundaryConditionType::Neumann;
	}

	//--------------//
	//     Misc     //
	//--------------//

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

		os << ": ";
		Shape()->Serialize(os);
	}

	virtual ~Face() {}
};