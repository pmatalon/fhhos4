#pragma once
#include "../FunctionalBasis/FunctionalBasis.h"
#include "HHOParameters.h"

template <int Dim>
class Poisson_HHO_Element;

template <int Dim>
class Poisson_HHO_Face : virtual public Face<Dim>
{
private:
	DenseMatrix _faceMassMatrix;
	DenseMatrix _invFaceMassMatrix;

	DenseMatrix _elem1_massCellFace;
	DenseMatrix _elem1_massReconstructFace;
	DenseMatrix _elem2_massCellFace;
	DenseMatrix _elem2_massReconstructFace;

	// Project a (k+1)-polynomial on the face
	DenseMatrix _elem1_projFromReconstruct;
	DenseMatrix _elem2_projFromReconstruct;

	// Project a k-polynomial on the face
	DenseMatrix _elem1_projFromCell;
	DenseMatrix _elem2_projFromCell;
public:
	HHOParameters<Dim>* HHO;

	Poisson_HHO_Face(BigNumber number, Element<Dim>* element1, Element<Dim>* element2) : Face<Dim>(number, element1, element2) {}

	void InitHHO(HHOParameters<Dim>* hho)
	{
		if (this->_faceMassMatrix.rows() > 0)
			return;

		this->HHO = hho;

		this->_faceMassMatrix = this->FaceMassMatrix(HHO->FaceBasis);
		this->_invFaceMassMatrix = this->_faceMassMatrix.inverse();

		if (this->Element1 != nullptr)
		{
			this->_elem1_massReconstructFace = this->MassMatrix(HHO->FaceBasis, this->Element1, HHO->ReconstructionBasis);
			this->_elem1_massCellFace = this->MassMatrix(HHO->FaceBasis, this->Element1, HHO->CellBasis);
			this->_elem1_projFromReconstruct = this->_invFaceMassMatrix * _elem1_massReconstructFace;
			this->_elem1_projFromCell = this->_invFaceMassMatrix * _elem1_massCellFace;
		}

		if (this->Element2 != nullptr)
		{
			this->_elem2_massReconstructFace = this->MassMatrix(HHO->FaceBasis, this->Element2, HHO->ReconstructionBasis);
			this->_elem2_massCellFace = this->MassMatrix(HHO->FaceBasis, this->Element2, HHO->CellBasis);
			this->_elem2_projFromReconstruct = this->_invFaceMassMatrix * _elem2_massReconstructFace;
			this->_elem2_projFromCell = this->_invFaceMassMatrix * _elem2_massCellFace;
		}
	}

	DenseMatrix FaceMassMatrix()
	{
		return this->_faceMassMatrix;
	}

	DenseMatrix InvFaceMassMatrix()
	{
		return this->_invFaceMassMatrix;
	}

	DenseMatrix FaceMassMatrix(FunctionalBasis<Dim-1>* basis)
	{
		return this->Shape()->FaceMassMatrix(basis);
	}

	DenseMatrix MassMatrix(FunctionalBasis<Dim - 1>* basis, Element<Dim>* element, FunctionalBasis<Dim>* cellBasis)
	{
		DenseMatrix M(basis->LocalFunctions.size(), cellBasis->LocalFunctions.size());
		for (BasisFunction<Dim - 1>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : cellBasis->LocalFunctions)
			{
				double term = this->ComputeMassTerm(phi1, element, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
			}
		}
		return M;
	}

	double ComputeMassTerm(BasisFunction<Dim - 1>* facePhi, Element<Dim>* element, BasisFunction<Dim>* reconstructPhi)
	{
		auto reconstructPhiOnFace = element->EvalPhiOnFace(this, reconstructPhi);

		RefFunction functionToIntegrate = [facePhi, reconstructPhiOnFace](RefPoint p) {
			return facePhi->Eval(p) * reconstructPhiOnFace(p);
		};

		int polynomialDegree = facePhi->GetDegree() + reconstructPhi->GetDegree();
		return this->Integral(functionToIntegrate, polynomialDegree);
	}

	DenseMatrix GetMassCellFace(Element<Dim>* element)
	{
		if (element == this->Element1)
			return _elem1_massCellFace;
		else if (element == this->Element2)
			return _elem2_massCellFace;
		assert(false);
	}

	DenseMatrix GetMassReconstructFace(Element<Dim>* element)
	{
		if (element == this->Element1)
			return _elem1_massReconstructFace;
		else if (element == this->Element2)
			return _elem2_massReconstructFace;
		assert(false);
	}

	DenseMatrix GetProjFromReconstruct(Element<Dim>* element)
	{
		if (element == this->Element1)
			return _elem1_projFromReconstruct;
		else if (element == this->Element2)
			return _elem2_projFromReconstruct;
		assert(false);
	}

	DenseMatrix GetProjFromCell(Element<Dim>* element)
	{
		if (element == this->Element1)
			return _elem1_projFromCell;
		else if (element == this->Element2)
			return _elem2_projFromCell;
		assert(false);
	}

	DenseMatrix GetProjFromCell(Element<Dim>* element, FunctionalBasis<Dim>* cellInterpolationBasis)
	{
		DenseMatrix massFaceCell = this->MassMatrix(HHO->FaceBasis, element, cellInterpolationBasis);
		DenseMatrix projFromCell = this->_invFaceMassMatrix * massFaceCell;
		return projFromCell;
	}

	double ProjectOnBasisFunction(BasisFunction<Dim - 1>* phi, DomFunction f)
	{
		RefFunction functionToIntegrate = [this, f, phi](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return f(domainPoint) * phi->Eval(refElementPoint);
		};

		return this->Integral(functionToIntegrate);
	}

	Vector ProjectOnBasis(FunctionalBasis<Dim - 1>* faceBasis, DomFunction f)
	{
		Vector projection(faceBasis->Size());
		for (BasisFunction<Dim - 1>* phi : faceBasis->LocalFunctions)
		{
			projection(phi->LocalNumber) = ProjectOnBasisFunction(phi, f);
		}
		return projection;
	}

	DenseMatrix ComputeCanonicalInjectionMatrixCoarseToFine(FunctionalBasis<Dim-1>* faceBasis)
	{
		DenseMatrix J(faceBasis->Size() * this->FinerFaces.size(), faceBasis->Size());

		for (auto f : this->FinerFaces)
		{
			Poisson_HHO_Face<Dim>* fineFace = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

			DenseMatrix fineCoarseMass(faceBasis->Size(), faceBasis->Size());
			for (BasisFunction<Dim-1>* finePhi : faceBasis->LocalFunctions)
			{
				for (BasisFunction<Dim-1>* coarsePhi : faceBasis->LocalFunctions)
				{
					RefFunction functionToIntegrate = [this, fineFace, finePhi, coarsePhi](RefPoint fineRefPoint) {
						DomPoint domPoint = fineFace->ConvertToDomain(fineRefPoint);
						RefPoint coarseRefPoint = this->ConvertToReference(domPoint);
						return finePhi->Eval(fineRefPoint)*coarsePhi->Eval(coarseRefPoint);
					};

					int polynomialDegree = finePhi->GetDegree() + coarsePhi->GetDegree();
					double integral = fineFace->Integral(functionToIntegrate, polynomialDegree);
					fineCoarseMass(finePhi->LocalNumber, coarsePhi->LocalNumber) = integral;
				}
			}

			DenseMatrix invFineMass = fineFace->InvFaceMassMatrix();

			J.block(this->LocalNumberOf(fineFace)*faceBasis->Size(), 0, faceBasis->Size(), faceBasis->Size()) = invFineMass * fineCoarseMass;
		}

		return J;
	}
};