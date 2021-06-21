#pragma once
#include "../FunctionalBasis/FunctionalBasis.h"
#include "HHOParameters.h"

template <int Dim>
class Diff_HHOElement;

template <int Dim>
class Diff_HHOFace : virtual public Face<Dim>
{
private:
	DenseMatrix _massMatrix;
	DenseMatrix _invMassMatrix;
public:
	HHOParameters<Dim>* HHO;

	Diff_HHOFace() {}

	Diff_HHOFace(BigNumber number, Element<Dim>* element1, Element<Dim>* element2) : Face<Dim>(number, element1, element2) {}

	void InitHHO(HHOParameters<Dim>* hho)
	{
		if (this->_massMatrix.rows() > 0)
			return;

		this->HHO = hho;

		//this->ComputeAndSaveQuadraturePoints(HHO->FaceBasis->GetDegree());
		//this->ComputeAndSaveQuadraturePoints();

		//this->_faceMassMatrix = this->FaceMassMatrix(HHO->FaceBasis);
		this->_massMatrix = this->Shape()->MassMatrix(HHO->FaceBasis);
		this->_invMassMatrix = this->_massMatrix.inverse();
	}

	DenseMatrix MassMatrix()
	{
		return this->_massMatrix;
	}

	DenseMatrix InvMassMatrix()
	{
		return this->_invMassMatrix;
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
		return ComputeTraceMatrix(element, HHO->ReconstructionBasis);
	}

	DenseMatrix TraceUsingCellBasis(Element<Dim>* element)
	{
		return ComputeTraceMatrix(element, HHO->CellBasis);
	}

	DenseMatrix Trace(Element<Dim>* element, FunctionalBasis<Dim>* cellInterpolationBasis)
	{
		/*if (cellInterpolationBasis == HHO->CellBasis)
			return TraceUsingCellBasis(element);
		else if (cellInterpolationBasis == HHO->ReconstructionBasis)
			return TraceUsingReconstructBasis(element);
		//else*/
			return ComputeTraceMatrix(element, cellInterpolationBasis);
		//assert(false);
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

			DenseMatrix invFineMass = fineFace->InvMassMatrix();

			J.block(this->LocalNumberOf(fineFace)*faceBasis->Size(), 0, faceBasis->Size(), faceBasis->Size()) = invFineMass * fineCoarseMass;
		}

		return J;
	}

	void DeleteUselessMatricesAfterAssembly()
	{
		Utils::Empty(_massMatrix);
	}

	void DeleteUselessMatricesAfterMultigridSetup()
	{
		Utils::Empty(_invMassMatrix);
		this->EmptySavedDomPoints();
	}

private:
	DenseMatrix ComputeTraceMatrix(Element<Dim>* e, FunctionalBasis<Dim>* cellBasis)
	{
		assert(_invMassMatrix.rows() > 0);
		DenseMatrix massCellFace = this->MassMatrix(HHO->FaceBasis, e, cellBasis);
		return _invMassMatrix * massCellFace;
	}
};