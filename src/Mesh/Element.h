#pragma once
#include <map>
#include <set>
#include <mutex>
#include "Vertex.h"
#include "../Geometry/PhysicalShape.h"
#include "../TestCases/Diffusion/DiffusionField.h"
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
	bool IsFullyEmbeddedInCoarseElement = true;
	vector<Face<Dim>*> FinerFacesRemoved;

	// Used for non-nested meshes. It contains at least FinerElements.
	map<Element<Dim>*, vector<PhysicalShape<Dim>*>> OverlappingFineElements;

	// Used during the non-nested mesh coarsening.
	// Coarse faces not to cross during element refinement for the L2 projection.
	set<Face<Dim>*> CoarseFacesNotToCross;

	// Used during mesh construction
	mutex Mutex;
	bool IsDeleted = false;
	bool HasReEntrantCorner = false;

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
	virtual void ExportToMatlab(string color = "r", bool writeNumber = false) const
	{
		Shape()->ExportToMatlab(color);
		if (writeNumber)
		{
			MatlabScript s;
			s.PlotText(Center(), to_string(this->Number));
		}
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

	virtual void RefineWithoutCoarseOverlap()
	{
		vector<PhysicalShape<Dim - 1>*> doNotCross;
		for (Face<Dim>* cf : CoarseFacesNotToCross)
		{
			if (!cf->IsDeleted)
				doNotCross.push_back(cf->Shape());
			else // should never happen
				Utils::Error("Coarse face " + to_string(cf->Number) + " deleted (fine elem id " + to_string(this->Id) + "). This should never happen.");
		}
		Shape()->RefineWithoutCoarseOverlap(doNotCross);

		CoarseFacesNotToCross = set<Face<Dim>*>(); // free memory
	}

	virtual void Refine(int nRefinements)
	{
		assert(false && "To be implemented");
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

	void AddOverlappingFineSubShape(Element<Dim>* fineElement, PhysicalShape<Dim>* subShape)
	{
		auto it = this->OverlappingFineElements.find(fineElement);
		if (it == this->OverlappingFineElements.end())
			this->OverlappingFineElements.insert({ fineElement , {subShape} });
		else
			it->second.push_back(subShape);
	}


	inline int LocalNumberOf(Element<Dim>* finerElement)
	{
		return this->_finerElementsLocalNumbering[finerElement];
	}
	inline int LocalNumberOfOverlapping(Element<Dim>* finerElement)
	{
		auto it = this->_overlappingFineElementsLocalNumbering.find(finerElement);
		if (it != this->_overlappingFineElementsLocalNumbering.end())
			return it->second;
		else
		{
			assert(false && "_overlappingFineElementsLocalNumbering probably not initialized");
			return -1;
		}
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

	bool HasAny(const vector<Face<Dim>*>& faces)
	{
		for (Face<Dim>* f : faces)
		{
			if (this->HasFace(f))
				return true;
		}
		return false;
	}

	bool HasAny(const vector<Vertex*>& vertices)
	{
		for (Vertex* v : vertices)
		{
			if (this->HasVertex(v, false))
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
		for (Vertex* v : this->Vertices())
		{
			if (other->HasVertex(v))
				return true;
		}
		return false;
	}

	bool IsIn(const vector<Element<Dim>*>& list)
	{
		auto it = find(list.begin(), list.end(), this);
		return it != list.end();
	}

	

public:
	virtual void RemoveIntersections(const vector<Vertex*>& verticesToRemove)
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

		set<Element<Dim>*> tested;

		set<Element<Dim>*> remaining;
		remaining.insert(faceNeigbours.begin(), faceNeigbours.end());
		while (!remaining.empty())
		{
			auto it = remaining.begin();
			Element<Dim>* candidate = *it;

			remaining.erase(it);

			if (candidate->HasOneVertexInCommonWith(this))
			{
				vertexNeighbours.insert(candidate);

				for (Element<Dim>* n : candidate->Neighbours(onlyInTheSamePhysicalPart))
				{
					if (n == this || tested.find(n) != tested.end())
						continue;
					remaining.insert(n);
				}
			}
			tested.insert(candidate);
		}

		return vector<Element<Dim>*>(vertexNeighbours.begin(), vertexNeighbours.end());
	}

	vector<Element<Dim>*> ThisAndVertexNeighbours(bool onlyInTheSamePhysicalPart = true)
	{
		vector<Element<Dim>*> list{ this };
		return Utils::Join(list, VertexNeighbours(onlyInTheSamePhysicalPart));
	}

	bool IsInSamePhysicalPartAs(Element<Dim>* other)
	{
		return (!this->PhysicalPart && !other->PhysicalPart) || this->PhysicalPart == other->PhysicalPart;
	}

	bool IsAtPhysicalPartBoundary()
	{
		for (Element<Dim>* n : this->Neighbours(false))
		{
			if (!this->IsInSamePhysicalPartAs(n))
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

	inline DenseMatrix MassMatrix(FunctionalBasis<Dim>* basis) const
	{
		return this->Shape()->MassMatrix(basis);
	}

	double IntegralKGradGrad(const Tensor<Dim>& K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return this->Shape()->IntegralKGradGrad(K, phi1, phi2);
	}

	//-----------------------------//
	//             HHO             //
	//-----------------------------//

	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis) const
	{
		return this->Shape()->CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	//--------------//
	//     Misc     //
	//--------------//

	// Constant diffusion tensor
	const Tensor<Dim>& DiffTensor() const
	{
		assert(PhysicalPart && "This element has no physical part.");
		assert(PhysicalPart->ConstantDiffTensor && "This element has no tensor.");
		return *PhysicalPart->ConstantDiffTensor;
	}

	// Constant isotropic diffusion coefficient (deprecated but still used in DG)
	double Kappa() const
	{
		return DiffTensor().LargestEigenValue;
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
		if (this->DiffTensor().LargestEigenValue > 1)
		{
			os << ", k=";
			os << this->DiffTensor().LargestEigenValue;
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
				if (this->Shape()->IsGeneralPolygon())
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