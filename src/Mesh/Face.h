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

	virtual PhysicalShape<Dim-1>* Shape() const = 0;

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
	virtual vector<Vertex*> Vertices() const
	{
		return Shape()->Vertices();
	}
	virtual void ExportToMatlab(string color = "r") const
	{
		return Shape()->ExportToMatlab(color);
	}

	bool HasVertex(Vertex* v, bool compareCoordinates = false)
	{
		return this->Shape()->HasVertex(v, compareCoordinates);
	}

	// Transformation to reference element
	virtual DomPoint ConvertToDomain(RefPoint refPoint) const
	{
		return Shape()->ConvertToDomain(refPoint);
	}
	virtual RefPoint ConvertToReference(DomPoint domainPoint) const
	{
		return Shape()->ConvertToReference(domainPoint);
	}

	//---------------------------------------//
	//   Relations to other faces/elements   //
	//---------------------------------------//

	bool HasSameVertices(Face<Dim>* other, bool compareCoordinates = false)
	{
		return this->Shape()->HasSameVertices(other->Shape(), compareCoordinates);
	}

	bool HasCommonVerticesWith(Face<Dim>* other, bool compareCoordinates = false)
	{
		for (Vertex* v : this->Vertices())
		{
			if (other->HasVertex(v, compareCoordinates))
				return true;
		}
		return false;
	}

	bool IntersectsWith(Face<Dim>* other)
	{
		assert(Dim == 2);

		// https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection

		DomPoint P1 = *this->Vertices()[0];
		double x1 = P1.X;
		double y1 = P1.Y;
		DomPoint P2 = *this->Vertices()[1];
		double x2 = P2.X;
		double y2 = P2.Y;

		DomPoint P3 = *other->Vertices()[0];
		double x3 = P3.X;
		double y3 = P3.Y;
		DomPoint P4 = *other->Vertices()[1];
		double x4 = P4.X;
		double y4 = P4.Y;

		double denom = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4);
		if (abs(denom) < Utils::NumericalZero)
			return false; // the two lines are parallel

		double t =  ((x1 - x3)*(y3 - y4) - (y1 - y3)*(x3 - x4)) / denom;
		double u = -((x1 - x2)*(y1 - y3) - (y1 - y2)*(x1 - x3)) / denom;

		double eps = Utils::Eps;
		return t + eps > 0 && t < 1 + eps  // t in [0,1]
			&& u + eps > 0 && u < 1 + eps; // u in [0,1]
	}

	bool IsIn(vector<Face<Dim>*> list)
	{
		for (Face<Dim>* f : list)
		{
			if (f == this)
				return true;
		}
		return false;
	}

	static int NumberOfFacesContainingVertex(vector<Face<Dim>*> faces, Vertex* v)
	{
		int nFacesContainingV = 0;
		for (Face<Dim>* f : faces)
		{
			if (f->HasVertex(v, true))
			{
				nFacesContainingV++;
				if (nFacesContainingV >= 2)
					break;
			}
		}
		return nFacesContainingV;
	}

	static bool IsInFaces(vector<Face<Dim>*> faces, Vertex* v)
	{
		bool vertexIsInFaces = false;
		for (Face<Dim>* f : faces)
		{
			if (f->HasVertex(v, true))
			{
				vertexIsInFaces = true;
				break;
			}
		}
		return vertexIsInFaces;
	}

	static bool IsInTwoFaces(vector<Face<Dim>*> faces, Vertex* v)
	{
		return NumberOfFacesContainingVertex(faces, v) >= 2;
	}

	virtual Face<Dim>* CreateSameGeometricFace(BigNumber number, Element<Dim>* element1) = 0;

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

	bool HasBeenCoarsened()
	{
		return this->CoarseFace->FinerFaces.size() > 1;
	}

	//------------------------//
	//        Integral        //
	//------------------------//

	virtual double Integral(RefFunction func) const
	{
		return Shape()->Integral(func);
	}
	virtual double Integral(RefFunction func, int polynomialDegree) const
	{
		return Shape()->Integral(func, polynomialDegree);
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