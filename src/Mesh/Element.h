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
public:
	BigNumber Id = 0;
	BigNumber Number;
	vector<Face<Dim>*> Faces;

	PhysicalGroup<Dim>* PhysicalPart = nullptr;

	// Intergrid links //
	vector<Element<Dim>*> FinerElements;
	Element<Dim>* CoarserElement = nullptr;
	vector<Face<Dim>*> FinerFacesRemoved;
	vector<Element<Dim>*> OverlappingFineElements; // Used for non-nested meshes. It contains at least FinerElements.
	bool IsFullyEmbeddedInCoarseElement = false;

	mutex Mutex;
	bool IsDeleted = false;

	Element(BigNumber number)
	{
		this->Number = number;
	}

	//-----------------------//
	//   Virtual functions   //
	//-----------------------//

	virtual PhysicalShape<Dim>* Shape() const
	{
		if (this->IsDeleted)
			Utils::FatalError("A method has been called on a deleted element!");
		else
			Utils::FatalError("The method Element.Shape() should be implemented in the subclass.");
		return nullptr;
	}

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
	virtual vector<Vertex*> Vertices() const
	{
		return Shape()->Vertices();
	}
	virtual bool Contains(DomPoint p) const
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
	virtual DomPoint ConvertToDomain(RefPoint refPoint) const
	{
		return Shape()->ConvertToDomain(refPoint);
	}
	virtual RefPoint ConvertToReference(DomPoint domainPoint) const
	{
		return Shape()->ConvertToReference(domainPoint);
	}
private:
	inline DimMatrix<Dim> InverseJacobianTranspose(RefPoint p) const
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
		for (Face<Dim>* f : this->Faces)
		{
			if (f == face)
				return true;
		}
		return false;
	}

	bool HasVertex(Vertex* v, bool compareCoordinates = false)
	{
		return this->Shape()->HasVertex(v, compareCoordinates);
	}

	bool HasSameVertices(Element<Dim>* other, bool compareCoordinates = false)
	{
		return this->Shape()->HasSameVertices(other->Shape(), compareCoordinates);
	}

	bool IsIn(vector<Element<Dim>*> list)
	{
		for (Element<Dim>* e : list)
		{
			if (e == this)
				return true;
		}
		return false;
	}

	// Replace faces with their agglomeration //
	void ReplaceFaces(vector<Face<Dim>*> faces, Face<Dim>* collapsedFace)
	{
		vector<Face<Dim>*> currentFaces(this->Faces);
		this->Faces.clear();
		for (Face<Dim>* f : currentFaces)
		{
			if (!f->IsIn(faces))
				this->Faces.push_back(f);
		}
		this->Faces.push_back(collapsedFace);

		RemoveIntersections(faces, collapsedFace);
	}

protected:
	virtual void RemoveIntersections(vector<Face<Dim>*> oldFaces, Face<Dim>* newFace)
	{
		assert(false);
	}
public:
	int LocalNumberOf(Element<Dim>* finerElement)
	{
		for (int i = 0; i < this->FinerElements.size(); i++)
		{
			if (this->FinerElements[i] == finerElement)
				return i;
		}
		assert(false);
	}

	int LocalNumberOfOverlapping(Element<Dim>* finerElement)
	{
		for (int i = 0; i < this->OverlappingFineElements.size(); i++)
		{
			if (this->OverlappingFineElements[i] == finerElement)
				return i;
		}
		assert(false);
	}

	inline int LocalNumberOf(Face<Dim>* face)
	{
		return this->_facesLocalNumbering[face];
	}

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

	vector<Element<Dim>*> Neighbours()
	{
		set<Element<Dim>*> neighbours;
		for (Face<Dim>* f : this->Faces)
		{
			if (f->IsDomainBoundary || f->IsDeleted)
				continue;
			Element<Dim>* neighbour = f->GetNeighbour(this);
			if (neighbour) // no neighbour might be affected yet
				neighbours.insert(neighbour);
		}

		return vector<Element<Dim>*>(neighbours.begin(), neighbours.end());
	}

	bool IsInSamePhysicalPartAs(Element<Dim>* other)
	{
		return (!this->PhysicalPart && !other->PhysicalPart) || this->PhysicalPart == other->PhysicalPart;
	}

	void CheckIfFullyEmbeddedInCoarseElement(CoarseningStrategy stgy)
	{
		if (Utils::BuildsNestedMeshHierarchy(stgy))
			this->IsFullyEmbeddedInCoarseElement = true;
		else
			this->IsFullyEmbeddedInCoarseElement == this->CoarserElement->FullyEmbeds(this);
	}

	void SetOverlappingFineElements(CoarseningStrategy stgy)
	{
		//if (!this->OverlappingFineElements.empty())
			//return;

		set<Element<Dim>*> overlapping;
		set<Element<Dim>*> tested;
		for (auto fe : this->FinerElements)
		{
			assert(fe->PhysicalPart == this->PhysicalPart);
			overlapping.insert(fe);
			tested.insert(fe);
		}

		if (!Utils::BuildsNestedMeshHierarchy(stgy))
		{
			for (Element<Dim>* neighbour : this->Neighbours())
			{
				for (auto fe : neighbour->FinerElements)
				{
					if (fe->PhysicalPart == this->PhysicalPart && !fe->IsFullyEmbeddedInCoarseElement && tested.find(fe) == tested.end() && fe->Overlaps(this))
						overlapping.insert(fe);
					tested.insert(fe);
					if (stgy != CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours)
					{
						for (Element<Dim>* feNeighbour : fe->Neighbours())
						{
							if (feNeighbour->PhysicalPart == this->PhysicalPart && !fe->IsFullyEmbeddedInCoarseElement && tested.find(feNeighbour) == tested.end() && feNeighbour->Overlaps(this))
								overlapping.insert(feNeighbour);
							tested.insert(feNeighbour);
						}
					}
				}
			}
		}

		this->OverlappingFineElements = vector<Element<Dim>*>(overlapping.begin(), overlapping.end());
	}

	bool FullyEmbeds(Element<Dim>* other)
	{
		for (Vertex* v : other->Vertices())
		{
			if (!this->Contains(*v))
				return false;
		}
		return true;
	}

	bool Overlaps(Element<Dim>* other)
	{
		for (DomPoint domPoint : this->QuadraturePoints())
		{
			if (other->Contains(domPoint))
				return true;
		}
		return false;
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
		RefFunction evalOnFace = [this, face, phi](RefPoint refPoint1D) {
			DomPoint domainPoint2D = face->ConvertToDomain(refPoint1D);
			RefPoint refPoint2D = this->ConvertToReference(domainPoint2D);
			return phi->Eval(refPoint2D);
		};
		return evalOnFace;
	}

	virtual function<DimVector<Dim>(RefPoint)> GradPhiOnFace(Face<Dim>* face, BasisFunction<Dim>* phi)
	{
		function<DimVector<Dim>(RefPoint)> gradOnFace = [this, face, phi](RefPoint refPoint1D) {
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
		RefFunction errorFunction = [this, exactSolution, approximate](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return pow(exactSolution(domainPoint) - approximate(refElementPoint), 2);
		};

		return Integral(errorFunction);
	}

	double EvalApproximateSolution(FunctionalBasis<Dim>* basis, const Vector &globalCoeffs, DomPoint evaluationPoint)
	{
		RefFunction localApproximate = basis->GetApproximateFunction(globalCoeffs, this->Number * basis->NumberOfLocalFunctionsInElement(this));
		return localApproximate(this->ConvertToReference(evaluationPoint));
	}

	double SourceTerm(BasisFunction<Dim>* phi, DomFunction f)
	{
		RefFunction sourceTimesBasisFunction = [this, f, phi](RefPoint refElementPoint) {
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
				for (Vertex* v : f->Shape()->Vertices())
				{
					DimVector<Dim> VC = Vect<Dim>(*v, C);
					assert(VC.dot(n) < 0);

					for (Vertex* v2 : f->Shape()->Vertices())
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