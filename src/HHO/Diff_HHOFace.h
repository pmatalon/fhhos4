#pragma once
#include "../FunctionalBasis/FunctionalBasis.h"
#include "HHOParameters.h"

template <int Dim>
class Diff_HHOElement;

template <int Dim>
class Diff_HHOFace : virtual public Face<Dim>
{
private:
	DenseMatrix _faceMassMatrix;
	DenseMatrix _invFaceMassMatrix;

	// Takes the trace of degree k of the reconstruction basis of element 1
	DenseMatrix _traceFromElem1UsingReconstructBasis;
	// Takes the trace of degree k of the reconstruction basis of element 2
	DenseMatrix _traceFromElem2UsingReconstructBasis;

	// Takes the trace of cell basis of element 1
	DenseMatrix _traceFromElem1UsingCellBasis;
	// Takes the trace of cell basis of element 2
	DenseMatrix _traceFromElem2UsingCellBasis;
public:
	HHOParameters<Dim>* HHO;

	Diff_HHOFace(BigNumber number, Element<Dim>* element1, Element<Dim>* element2) : Face<Dim>(number, element1, element2) {}

	void InitHHO(HHOParameters<Dim>* hho)
	{
		if (this->_faceMassMatrix.rows() > 0)
			return;

		this->HHO = hho;

		//this->ComputeAndSaveQuadraturePoints(HHO->FaceBasis->GetDegree());
		this->ComputeAndSaveQuadraturePoints();

		this->_faceMassMatrix = this->FaceMassMatrix(HHO->FaceBasis);
		this->_invFaceMassMatrix = this->_faceMassMatrix.inverse();

		if (this->Element1 != nullptr)
		{
			this->_traceFromElem1UsingCellBasis        = ComputeTraceMatrix(this->Element1, HHO->CellBasis);
			this->_traceFromElem1UsingReconstructBasis = ComputeTraceMatrix(this->Element1, HHO->ReconstructionBasis);
		}

		if (this->Element2 != nullptr)
		{
			this->_traceFromElem2UsingCellBasis        = ComputeTraceMatrix(this->Element2, HHO->CellBasis);
			this->_traceFromElem2UsingReconstructBasis = ComputeTraceMatrix(this->Element2, HHO->ReconstructionBasis);
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

		RefFunction functionToIntegrate = [facePhi, reconstructPhiOnFace](const RefPoint& p) {
			return facePhi->Eval(p) * reconstructPhiOnFace(p);
		};

		int polynomialDegree = facePhi->GetDegree() + reconstructPhi->GetDegree();
		return this->Integral(functionToIntegrate, polynomialDegree);
	}

	DenseMatrix TraceUsingReconstructBasis(Element<Dim>* element)
	{
		assert(_traceFromElem1UsingReconstructBasis.rows() > 0);
		if (element == this->Element1)
			return _traceFromElem1UsingReconstructBasis;
		else if (element == this->Element2)
			return _traceFromElem2UsingReconstructBasis;
		else
			return ComputeTraceMatrix(element, HHO->ReconstructionBasis);
	}

	DenseMatrix TraceUsingCellBasis(Element<Dim>* element)
	{
		assert(_traceFromElem1UsingCellBasis.rows() > 0);
		if (element == this->Element1)
			return _traceFromElem1UsingCellBasis;
		else if (element == this->Element2)
			return _traceFromElem2UsingCellBasis;
		else
			return ComputeTraceMatrix(element, HHO->CellBasis);
	}

	DenseMatrix Trace(Element<Dim>* element, FunctionalBasis<Dim>* cellInterpolationBasis)
	{
		if (cellInterpolationBasis == HHO->CellBasis)
			return TraceUsingCellBasis(element);
		else if (cellInterpolationBasis == HHO->ReconstructionBasis)
			return TraceUsingReconstructBasis(element);
		//else
			//return ComputeTraceMatrix(element, cellInterpolationBasis);
		assert(false);
	}

	double ProjectOnBasisFunction(BasisFunction<Dim - 1>* phi, DomFunction f)
	{
		RefFunction functionToIntegrate = [this, f, phi](const RefPoint& refElementPoint) {
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
			Diff_HHOFace<Dim>* fineFace = dynamic_cast<Diff_HHOFace<Dim>*>(f);

			DenseMatrix fineCoarseMass(faceBasis->Size(), faceBasis->Size());
			for (BasisFunction<Dim-1>* finePhi : faceBasis->LocalFunctions)
			{
				for (BasisFunction<Dim-1>* coarsePhi : faceBasis->LocalFunctions)
				{
					RefFunction functionToIntegrate = [this, fineFace, finePhi, coarsePhi](const RefPoint& fineRefPoint) {
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

	void DeleteUselessMatricesAfterAssembly()
	{
		Utils::Empty(_faceMassMatrix);
	}

	void DeleteUselessMatricesAfterMultigridSetup()
	{
		Utils::Empty(_invFaceMassMatrix);
		Utils::Empty(_traceFromElem1UsingReconstructBasis);
		Utils::Empty(_traceFromElem2UsingReconstructBasis);
		Utils::Empty(_traceFromElem1UsingCellBasis);
		Utils::Empty(_traceFromElem2UsingCellBasis);
		this->EmptySavedDomPoints();
	}

private:
	DenseMatrix ComputeTraceMatrix(Element<Dim>* e, FunctionalBasis<Dim>* cellBasis)
	{
		assert(_invFaceMassMatrix.rows() > 0);
		DenseMatrix massCellFace = this->MassMatrix(HHO->FaceBasis, e, cellBasis);
		return _invFaceMassMatrix * massCellFace;
	}
};