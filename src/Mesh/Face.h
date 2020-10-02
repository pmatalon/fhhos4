#pragma once
#include "Element.h"
#include "../Problem/BoundaryConditions.h"
using namespace std;

template <int Dim>
class Face
{
public:
	BigNumber Number;
	bool IsDomainBoundary;
	BoundaryGroup* BoundaryPart = nullptr;
	Element<Dim>* Element1;
	Element<Dim>* Element2;

	bool IsRemovedOnCoarserGrid = false; 
	vector<Face<Dim>*> FinerFaces;
	Face<Dim>* CoarseFace = nullptr;

	mutex Mutex;
	bool IsDeleted = false;
public:
	Face() { assert(false); }
	Face(BigNumber number, Element<Dim>* element1, Element<Dim>* element2)
	{
		this->Number = number;
		this->Element1 = element1;
		this->Element2 = element2;
		this->IsDomainBoundary = (element2 == nullptr);
	}
	Face(BigNumber number, Element<Dim>* element1)
		:Face(number, element1, nullptr)
	{
		this->IsDomainBoundary = true;
	}

	//-----------------------//
	//   Virtual functions   //
	//-----------------------//

	virtual const PhysicalShape<Dim - 1>* Shape() const = 0;
	/*{
		if (this->IsDeleted)
			Utils::FatalError("A method has been called on a deleted face!");
		else
			Utils::FatalError("The method Face.Shape() should be implemented in the subclass.");
		return nullptr;
	}*/
	virtual PhysicalShape<Dim - 1>* Shape() = 0;
	/*{
		if (this->IsDeleted)
			Utils::FatalError("A method has been called on a deleted face!");
		else
			Utils::FatalError("The method Face.Shape() should be implemented in the subclass.");
		return nullptr;
	}*/

	virtual vector<Vertex*> Vertices() const = 0;

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
	virtual bool Contains(const DomPoint& p) const
	{
		return Shape()->Contains(p);
	}
	virtual void ExportToMatlab(string color = "r") const
	{
		return Shape()->ExportToMatlab(color);
	}

	// Transformation to reference element
	virtual DomPoint ConvertToDomain(const RefPoint& refPoint) const
	{
		return Shape()->ConvertToDomain(refPoint, true);
	}
	virtual RefPoint ConvertToReference(const DomPoint& domainPoint) const
	{
		return Shape()->ConvertToReference(domainPoint);
	}

	//----------------------------------------//
	// Correspondance RefPoint/DomPoint saved //
	//----------------------------------------//

	DomPoint ConvertToDomainAndSaveResult(const RefPoint& refPoint)
	{
		return Shape()->ConvertToDomainAndSaveResult(refPoint);
	}
	void ComputeAndSaveQuadraturePoints(int polynomialDegree)
	{
		Shape()->ComputeAndSaveQuadraturePoints(polynomialDegree);
	}
	void ComputeAndSaveQuadraturePoints()
	{
		Shape()->ComputeAndSaveQuadraturePoints();
	}
	inline void EmptySavedDomPoints()
	{
		Shape()->EmptySavedDomPoints();
	}

	//---------------------------------------//
	//   Relations to other faces/elements   //
	//---------------------------------------//

	bool HasVertex(Vertex* v, bool compareCoordinates = false)
	{
		auto vertices = this->Vertices();
		auto it = find_if(vertices.begin(), vertices.end(), [v, compareCoordinates](Vertex* v2) { return v == v2 || (compareCoordinates && *v == *v2); });
		return it != vertices.end();
	}

	bool HasSameVertices(Face<Dim>* other, bool compareCoordinates = false)
	{
		if (this->Vertices().size() != other->Vertices().size())
			return false;

		for (Vertex* v : this->Vertices())
		{
			if (!other->HasVertex(v, compareCoordinates))
				return false;
		}
		return true;
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

	bool IsIn(const vector<Face<Dim>*>& list)
	{
		auto it = find(list.begin(), list.end(), this);
		return it != list.end();
	}

	static int NumberOfFacesContainingVertex(const vector<Face<Dim>*>& faces, Vertex* v)
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

	static bool IsInFaces(const vector<Face<Dim>*>& faces, Vertex* v)
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

	static bool IsInTwoFaces(const vector<Face<Dim>*>& faces, Vertex* v)
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
		return nullptr;
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

	Face<Dim>* ClosestFaceAmongst(const vector<Face<Dim>*>& list, bool isDomainBoundary)
	{
		double smallestDistance = -1;
		Face<Dim>* closestFace = nullptr;
		for (Face<Dim>* f : list)
		{
			if (f->IsDomainBoundary == isDomainBoundary)
			{
				double distance = Vect<Dim>(this->Center(), f->Center()).norm();
				if (!closestFace || distance < smallestDistance)
				{
					smallestDistance = distance;
					closestFace = f;
				}
			}
		}
		return closestFace;
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

	inline bool HasDirichletBC()
	{
		return IsDomainBoundary && (!BoundaryPart || BoundaryPart->Condition == BoundaryConditionType::Dirichlet);
	}

	inline bool HasNeumannBC()
	{
		return IsDomainBoundary && BoundaryPart && BoundaryPart->Condition == BoundaryConditionType::Neumann;
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

	virtual ~Face()
	{
		this->IsDeleted = true;
	}
};