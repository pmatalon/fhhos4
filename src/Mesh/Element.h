#pragma once
#include <vector>
#include <map>
#include "Vertex.h"
#include "../Utils/Utils.h"
#include "../Problem/SourceFunction.h"
#include "../Problem/DiffusionPartition.h"
#include "../FunctionalBasis/FunctionalBasis.h"

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

	// Geometric information
	virtual double Diameter() = 0;
	virtual double Measure() = 0;
	virtual DomPoint Center() = 0;
	virtual DimVector<Dim> OuterNormalVector(Face<Dim>* face) = 0;

	// Transformation to reference element
	virtual DomPoint ConvertToDomain(RefPoint refPoint) const = 0;
	virtual RefPoint ConvertToReference(DomPoint domainPoint) const = 0;
	virtual DimMatrix<Dim> InverseJacobianTranspose() const = 0;

	// Integral
	virtual double Integral(BasisFunction<Dim>* phi) const = 0;
	virtual double Integral(RefFunction func) const = 0;
	virtual double Integral(RefFunction func, int polynomialDegree) const = 0;

	//---------------------------//
	//   Geometric information   //
	//---------------------------//

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

	virtual double Integral(DomFunction globalFunction) const
	{
		RefFunction refFunction = [this, globalFunction](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return globalFunction(domainPoint);
		};

		return Integral(refFunction);
	}

	virtual double Integral(DomFunction globalFunction, int polynomialDegree) const
	{
		RefFunction refFunction = [this, globalFunction](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return globalFunction(domainPoint);
		};

		return Integral(refFunction, polynomialDegree);
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
			DimMatrix<Dim> invJ = this->InverseJacobianTranspose();
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
	}

	friend ostream& operator<<(ostream& os, const Element<Dim>& s)
	{
		s.Serialize(os);
		return os;
	}

	virtual ~Element() {}
};