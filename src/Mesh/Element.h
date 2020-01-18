#pragma once
#include <vector>
#include <map>
#include "Vertex.h"
#include "GeometricShapeWithReferenceShape.h"
#include "../Problem/SourceFunction.h"
#include "../Problem/DiffusionPartition.h"

template <int Dim>
class Face;


template <int Dim>
class Element 
{
private:
	std::map<Face<Dim>*, int> _facesLocalNumbering;
public:
	BigNumber Number;
	std::vector<Face<Dim>*> Faces;

	double Kappa = 1; // constant diffusion coefficient (used in DG)
	Tensor<Dim>* DiffTensor = nullptr;

	std::vector<Element<Dim>*> FinerElements;
	Element<Dim>* CoarserElement;
	std::vector<Face<Dim>*> FinerFacesRemoved;

	Element(BigNumber number)
	{
		this->Number = number;
	}

	//-----------------------//
	//   Virtual functions   //
	//-----------------------//

	virtual GeometricShapeWithReferenceShape<Dim>* Shape() const = 0;

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

	Element<Dim>* ElementOnTheOtherSideOf(Face<Dim>* face)
	{
		return face->GetNeighbour(this);
	}
	
	Face<Dim>* InterfaceWith(Element<Dim>* other)
	{
		for (Face<Dim>* f1 : this->Faces)
		{
			for (Face<Dim>* f2 : other->Faces)
			{
				if (f1 == f2)
					return f1;
			}
		}
		return nullptr;
	}

	int LocalNumberOf(Element<Dim>* finerElement)
	{
		for (int i = 0; i < this->FinerElements.size(); i++)
		{
			if (this->FinerElements[i] == finerElement)
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

	//-------------------//
	//     Integrals     //
	//-------------------//

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

	double L2ErrorPow2(RefFunction approximate, DomFunction exactSolution, int degreeUsed) const
	{
		RefFunction errorFunction = [this, exactSolution, approximate](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return pow(exactSolution(domainPoint) - approximate(refElementPoint), 2);
		};

		return Integral(errorFunction, degreeUsed);
	}

	double SourceTerm(BasisFunction<Dim>* phi, SourceFunction* f)
	{
		RefFunction sourceTimesBasisFunction = [this, f, phi](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return f->Eval(domainPoint) * phi->Eval(refElementPoint);
		};

		return Integral(sourceTimesBasisFunction);
	}

	//--------------//
	//     Misc     //
	//--------------//

	// For DG
	void SetDiffusionCoefficient(DiffusionPartition<Dim>* diffusionPartition)
	{
		this->Kappa = diffusionPartition->Coefficient(Center());
	}
	void SetDiffusionTensor(DiffusionPartition<Dim>* diffusionPartition)
	{
		this->DiffTensor = diffusionPartition->DiffTensor(Center());
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
	}

	friend ostream& operator<<(ostream& os, const Element<Dim>& s)
	{
		s.Serialize(os);
		return os;
	}

	virtual ~Element() {}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	virtual void UnitTests() const
	{
		Shape()->UnitTests();

		// The following tests are valid if the element is convex
		DomPoint C = this->Center();
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
};