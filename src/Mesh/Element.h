#pragma once
#include <map>
#include <set>
#include <mutex>
#include "Vertex.h"
#include "../Geometry/PhysicalShape.h"
#include "../Problem/DiffusionField.h"
#include "PhysicalGroup.h"
#include "../Utils/MatlabScript.h"
using namespace std;

template <int Dim>
class Face;


template <int Dim>
class Element 
{
private:
	map<Face<Dim>*, int> _facesLocalNumbering;
	map<Element<Dim>*, int> _finerElementsLocalNumbering;
	map<Element<Dim>*, int> _overlappingFineElementsLocalNumbering;
public:
	BigNumber Id = 0;
	BigNumber Number;
	vector<Face<Dim>*> Faces;

	PhysicalGroup<Dim>* PhysicalPart = nullptr;

	// Intergrid links //
	vector<Element<Dim>*> FinerElements;
	Element<Dim>* CoarserElement = nullptr;
	vector<Face<Dim>*> FinerFacesRemoved;
	// Used for non-nested meshes. It contains at least FinerElements.
	map<Element<Dim>*, vector<PhysicalShape<Dim>*>> OverlappingFineElements;

	// Used during mesh construction
	mutex Mutex;
	bool IsDeleted = false;

	//-----------------------//
	//      Constructors     //
	//-----------------------//

	Element() {}

	// Copy constructor
	Element(const Element<Dim>& e)
	{
		// The default copy constructor doesn't work because the mutex is not a copyable object
		Id = e.Id;
		Number = e.Number;
		PhysicalPart = e.PhysicalPart;
	}

	// Assignment operator
	Element<Dim>& operator=(const Element<Dim>& e)
	{
		Id = e.Id;
		Number = e.Number;
		PhysicalPart = e.PhysicalPart;
		return *this;
	}

	Element(BigNumber number)
	{
		this->Number = number;
	}

	//-----------------------//
	//   Virtual functions   //
	//-----------------------//

	virtual PhysicalShape<Dim>* Shape()
	{
		if (this->IsDeleted)
			Utils::FatalError("A method has been called on a deleted element!");
		else
			Utils::FatalError("The method Element.Shape() should be implemented in the subclass.");
		return nullptr;
	}

	virtual const PhysicalShape<Dim>* Shape() const
	{
		if (this->IsDeleted)
			Utils::FatalError("A method has been called on a deleted element!");
		else
			Utils::FatalError("The method Element.Shape() should be implemented in the subclass.");
		return nullptr;
	}

	virtual vector<Vertex*> Vertices() const = 0;

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
	virtual DomPoint InteriorPoint() const
	{
		return Shape()->InteriorPoint();
	}
	virtual bool IsConvex() const
	{
		return Shape()->IsConvex();
	}
	virtual double InRadius() const
	{
		return Shape()->InRadius();
	}
	virtual double Regularity() const
	{
		return Shape()->Regularity();
	}
	virtual bool Contains(const DomPoint& p) const
	{
		return Shape()->Contains(p);
	}
	virtual bool ConvexHullEmbeds(Element<Dim>* other) const
	{
		return Shape()->ConvexHullEmbeds(other->Shape());
	}
	virtual void ExportToMatlab(string color = "r") const
	{
		return Shape()->ExportToMatlab(color);
	}
	virtual DimVector<Dim> OuterNormalVector(Face<Dim>* face) const = 0;

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

	DomPoint ConvertToDomainAndSaveResult(const RefPoint& refPoint, bool lock = false)
	{
		return Shape()->ConvertToDomainAndSaveResult(refPoint, lock);
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
private:
	inline DimMatrix<Dim> InverseJacobianTranspose(const RefPoint& p) const
	{
		return Shape()->InverseJacobianTranspose(p);
	}

public:
	void AddFace(Face<Dim>* face)
	{
		this->Faces.push_back(face);

		int faceLocalNumber = static_cast<int>(this->_facesLocalNumbering.size());
		this->_facesLocalNumbering.insert(std::pair<Face<Dim>*, int>(face, faceLocalNumber));
	}

	void RemoveFace(Face<Dim>* f)
	{
		int i = 0;
		for (i = 0; i < this->Faces.size(); i++)
		{
			if (this->Faces[i] == f)
				break;
		}
		this->Faces.erase(this->Faces.begin() + i);
	}

	void InitFaceLocalNumbering()
	{
		this->_facesLocalNumbering.clear();
		int faceLocalNumber = 0;
		for (Face<Dim>* f : this->Faces)
			this->_facesLocalNumbering.insert({ f, faceLocalNumber++ });
	}
	void InitFinerElementsLocalNumbering()
	{
		this->_finerElementsLocalNumbering.clear();
		int localNumber = 0;
		for (Element<Dim>* e : this->FinerElements)
			this->_finerElementsLocalNumbering.insert({ e, localNumber++ });
	}
	void InitOverlappingElementsLocalNumbering()
	{
		this->_overlappingFineElementsLocalNumbering.clear();
		int localNumber = 0;
		for (auto it = this->OverlappingFineElements.begin(); it != this->OverlappingFineElements.end(); it++)
		{
			Element<Dim>* e = it->first;
			this->_overlappingFineElementsLocalNumbering.insert({ e, localNumber++ });
		}
	}


	inline int LocalNumberOf(Element<Dim>* finerElement)
	{
		return this->_finerElementsLocalNumbering[finerElement];
	}
	inline int LocalNumberOfOverlapping(Element<Dim>* finerElement)
	{
		return this->_overlappingFineElementsLocalNumbering[finerElement];
	}
	inline int LocalNumberOf(Face<Dim>* face)
	{
		return this->_facesLocalNumbering[face];
	}

	Element<Dim>* ElementOnTheOtherSideOf(Face<Dim>* face)
	{
		return face->GetNeighbour(this);
	}
	
	Face<Dim>* CommonFaceWith(Element<Dim>* other)
	{
		for (Face<Dim>* f1 : this->Faces)
		{
			if (f1->IsDomainBoundary)
				continue;
			for (Face<Dim>* f2 : other->Faces)
			{
				if (f1 == f2)
					return f1;
			}
		}
		return nullptr;
	}

	vector<Face<Dim>*> InterfaceWith(Element<Dim>* other)
	{
		vector<Face<Dim>*> faces;
		for (Face<Dim>* f : this->Faces)
		{
			if (f->IsDomainBoundary)
				continue;
			if (find(other->Faces.begin(), other->Faces.end(), f) != other->Faces.end())
				faces.push_back(f);
		}
		return faces;
	}

	vector<Face<Dim>*> NonCommonFacesWith(Element<Dim>* other)
	{
		return Utils::SymmetricDifference<Face<Dim>*>(this->Faces, other->Faces);
	}

	vector<Face<Dim>*> BoundaryFaces()
	{
		vector<Face<Dim>*> faces;
		for (Face<Dim>* f : this->Faces)
		{
			if (f->IsDomainBoundary)
				faces.push_back(f);
		}
		return faces;
	}

	bool HasFace(Face<Dim>* face)
	{
		return find(this->Faces.begin(), this->Faces.end(), face) != this->Faces.end();
	}

	bool HasAny(vector<Face<Dim>*> faces)
	{
		for (Face<Dim>* f : faces)
		{
			if (this->HasFace(f))
				return true;
		}
		return false;
	}

	bool HasVertex(Vertex* v, bool compareCoordinates = false)
	{
		auto vertices = this->Vertices();
		auto it = find_if(vertices.begin(), vertices.end(), [v, compareCoordinates](Vertex* v2) { return v == v2 || (compareCoordinates && *v == *v2); });
		return it != vertices.end();
	}

	bool HasSameVertices(Element<Dim>* other, bool compareCoordinates = false)
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

	bool HasOneVertexInCommonWith(Element<Dim>* other)
	{
		return this->Shape()->HasOneVertexInCommonWith(other->Shape());
	}

	bool IsIn(const vector<Element<Dim>*>& list)
	{
		auto it = find(list.begin(), list.end(), this);
		return it != list.end();
	}

	// Replace faces with their agglomeration //
	void ReplaceFaces(const vector<Face<Dim>*>& faces, Face<Dim>* collapsedFace)
	{
		vector<Face<Dim>*> currentFaces = this->Faces;
		this->Faces.clear();
		for (Face<Dim>* f : currentFaces)
		{
			if (!f->IsIn(faces))
				this->Faces.push_back(f);
		}
		this->Faces.push_back(collapsedFace);

		double oldMeasure = this->Measure();
		RemoveIntersections(faces, collapsedFace);
		double newMeasure = this->Measure();

		// Some fine elements associated to this coarse one might now be more overlapping the other coarse elements.
		// So we change the association.
		if (!collapsedFace->IsDomainBoundary && newMeasure < oldMeasure)
		{
			Element<Dim>* neighbour = collapsedFace->GetNeighbour(this);
			auto it = this->FinerElements.begin();
			while (it != this->FinerElements.end())
			{
				Element<Dim>* fe = *it;
				if (fe->HasAny(collapsedFace->FinerFaces))
				{
					DomPoint p = fe->InteriorPoint();
					if (!this->Contains(p) && !collapsedFace->Contains(p) && neighbour->Contains(p))
					{
						/*cout << "%---------- this " << endl;
						this->ExportToMatlab("b");
						cout << "% neighbour " << endl;
						neighbour->ExportToMatlab("m");
						cout << "% transfered finer element (b --> m) " << endl;
						fe->ExportToMatlab("r");
						MatlabScript s;
						s.PlotPoint(p, "rx");*/

						// Transfer to the neighbour
						fe->CoarserElement = neighbour;
						neighbour->FinerElements.push_back(fe);
						it = this->FinerElements.erase(it);
						continue;
					}
				}
				it++;
			}
		}

	}

protected:
	virtual void RemoveIntersections(const vector<Face<Dim>*>& oldFaces, Face<Dim>* newFace)
	{
		assert(false);
	}
public:

	bool IsOnBoundary()
	{
		for (Face<Dim>* f : Faces)
		{
			if (f->IsDomainBoundary)
				return true;
		}
		return false;
	}

	double FrontierMeasure() const
	{
		double measure = 0;
		for (auto f : this->Faces)
			measure += f->Measure();
		return measure;
	}

	vector<Element<Dim>*> Neighbours(bool onlyInTheSamePhysicalPart = false)
	{
		set<Element<Dim>*> neighbours;
		for (Face<Dim>* f : this->Faces)
		{
			if (f->IsDomainBoundary || f->IsDeleted)
				continue;
			Element<Dim>* neighbour = f->GetNeighbour(this);
			if (neighbour && (!onlyInTheSamePhysicalPart || this->IsInSamePhysicalPartAs(neighbour))) // no neighbour might be affected yet
				neighbours.insert(neighbour);
		}

		return vector<Element<Dim>*>(neighbours.begin(), neighbours.end());
	}

	vector<Element<Dim>*> NeighboursInSamePhysicalPart()
	{
		return Neighbours(true);
	}

	vector<Element<Dim>*> VertexNeighbours(bool onlyInTheSamePhysicalPart = true)
	{
		set<Element<Dim>*> vertexNeighbours;
		vector<Element<Dim>*> faceNeigbours = this->Neighbours(onlyInTheSamePhysicalPart);
		vertexNeighbours.insert(faceNeigbours.begin(), faceNeigbours.end());

		for (Element<Dim>* n : faceNeigbours)
		{
			for (Element<Dim>* n2 : n->Neighbours(onlyInTheSamePhysicalPart))
			{
				if (n2 == this)
					continue;

				for (Vertex* v : n2->Vertices())
				{
					if (this->HasVertex(v))
						vertexNeighbours.insert(n2);
				}
			}
		}

		return vector<Element<Dim>*>(vertexNeighbours.begin(), vertexNeighbours.end());
	}

	bool IsInSamePhysicalPartAs(Element<Dim>* other)
	{
		return (!this->PhysicalPart && !other->PhysicalPart) || this->PhysicalPart == other->PhysicalPart;
	}
	
	//-------------------//
	//     Integrals     //
	//-------------------//

	vector<DomPoint> QuadraturePoints() const
	{
		return Shape()->QuadraturePoints();
	}

	virtual double Integral(BasisFunction<Dim>* phi) const
	{
		return Shape()->Integral(phi);
	}
	virtual double Integral(RefFunction func) const
	{
		return Shape()->Integral(func);
	}
	virtual double Integral(RefFunction func, int polynomialDegree) const
	{
		return Shape()->Integral(func, polynomialDegree);
	}
	virtual double Integral(DomFunction globalFunction) const
	{
		return Shape()->Integral(globalFunction);
	}
	virtual double Integral(DomFunction globalFunction, int polynomialDegree) const
	{
		return Shape()->Integral(globalFunction, polynomialDegree);
	}

	virtual RefFunction EvalPhiOnFace(Face<Dim>* face, BasisFunction<Dim>* phi)
	{
		RefFunction evalOnFace = [this, face, phi](const RefPoint& refPoint1D) {
			DomPoint domainPoint2D = face->ConvertToDomain(refPoint1D);
			RefPoint refPoint2D = this->ConvertToReference(domainPoint2D);
			return phi->Eval(refPoint2D);
		};
		return evalOnFace;
	}

	virtual function<DimVector<Dim>(RefPoint)> GradPhiOnFace(Face<Dim>* face, BasisFunction<Dim>* phi)
	{
		function<DimVector<Dim>(RefPoint)> gradOnFace = [this, face, phi](const RefPoint& refPoint1D) {
			DomPoint domainPoint2D = face->ConvertToDomain(refPoint1D);
			RefPoint refPoint2D = this->ConvertToReference(domainPoint2D);
			DimVector<Dim> gradPhi = phi->Grad(refPoint2D);
			DimMatrix<Dim> invJ = this->InverseJacobianTranspose(refPoint2D);
			DimVector<Dim> result = invJ * gradPhi;
			return result;
		};
		return gradOnFace;
	}

	double L2ErrorPow2(RefFunction approximate, DomFunction exactSolution) const
	{
		RefFunction errorFunction = [this, exactSolution, approximate](const RefPoint& refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return pow(exactSolution(domainPoint) - approximate(refElementPoint), 2);
		};

		return Integral(errorFunction);
	}

	double EvalApproximateSolution(FunctionalBasis<Dim>* basis, const Vector &globalCoeffs, const DomPoint& evaluationPoint)
	{
		RefFunction localApproximate = basis->GetApproximateFunction(globalCoeffs, this->Number * basis->NumberOfLocalFunctionsInElement(this));
		return localApproximate(this->ConvertToReference(evaluationPoint));
	}

	double SourceTerm(BasisFunction<Dim>* phi, DomFunction f)
	{
		RefFunction sourceTimesBasisFunction = [this, f, phi](const RefPoint& refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return f(domainPoint) * phi->Eval(refElementPoint);
		};

		return Integral(sourceTimesBasisFunction);
	}

	//--------------//
	//     Misc     //
	//--------------//

	// Constant diffusion tensor
	Tensor<Dim>* DiffTensor() const
	{
		assert(PhysicalPart && "This element has no physical part.");
		return PhysicalPart->ConstantDiffTensor;
	}

	// Constant isotropic diffusion coefficient (deprecated but still used in DG)
	double Kappa() const
	{
		return DiffTensor()->LargestEigenValue;
	}

	virtual void Serialize(ostream& os) const
	{
		os << "Element " << this->Number << ": ";
		if (Dim == 1)
			os << "interfaces ";
		else if (Dim == 2)
			os << "edges ";
		else if (Dim == 3)
			os << "faces ";
		for (size_t i = 0; i < this->Faces.size(); ++i)
		{
			os << this->Faces[i]->Number;
			if (i < this->Faces.size() - 1)
				os << ", ";
		}
		os << ", ";
		Shape()->Serialize(os);
		if (this->DiffTensor() && this->DiffTensor()->LargestEigenValue > 1)
		{
			os << ", k=";
			os << this->DiffTensor()->LargestEigenValue;
		}
	}

	friend ostream& operator<<(ostream& os, const Element<Dim>& s)
	{
		s.Serialize(os);
		return os;
	}

	virtual ~Element()
	{ 
		this->IsDeleted = true;
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	virtual void UnitTests() const
	{
		Shape()->UnitTests();

		assert(this->PhysicalPart);

		// The following tests are valid if the element is convex
		if (this->IsConvex())
		{
			DomPoint C = this->Center();
			bool containsCenter = this->Contains(C); // comment this test if Contains() is not implemented
			if (!containsCenter)
			{
				if (this->Shape()->IsMadeOfSubShapes())
				{
					cout << "Analysis of the problem:" << endl;
					cout << "Matlab script to plot the subshapes of element " << this->Number << ":" << endl;
					this->Shape()->ExportSubShapesToMatlab();
					cout << endl;
				}
			}
			assert(containsCenter && "The element does not contain its centroid.");


			for (Face<Dim>* f : this->Faces)
			{
				DimVector<Dim> n = this->OuterNormalVector(f);
				for (Vertex* v : f->Vertices())
				{
					DimVector<Dim> VC = Vect<Dim>(*v, C);
					assert(VC.dot(n) < 0);

					for (Vertex* v2 : f->Vertices())
					{
						if (v2 != v)
						{
							DimVector<Dim> V1V2 = Vect<Dim>(v, v2);
							assert(abs(V1V2.dot(n)) < 1e-15);
						}
					}
				}
			}
		}

	}
};