#pragma once
#include "../FunctionalBasis/OrthonormalBasis.h"
#include "HHOParameters.h"

template <int Dim>
class Diff_HHOElement;

template <int Dim>
class Diff_HHOFace
{
private:
	DenseMatrix _massMatrix;
	DenseMatrix _invMassMatrix;
public:
	Face<Dim>* MeshFace = nullptr;
	HHOParameters<Dim>* HHO;
	FunctionalBasis<Dim - 1>* Basis;

	Diff_HHOFace() {}

	//--------------//
	// Face wrapper //
	//--------------//

public:
	inline BigNumber Number() const
	{
		return this->MeshFace->Number;
	}
	inline double Diameter() const
	{
		return this->MeshFace->Diameter();
	}
	inline bool HasDirichletBC()
	{
		return this->MeshFace->HasDirichletBC();
	}
	inline bool HasNeumannBC()
	{
		return this->MeshFace->HasNeumannBC();
	}
private:
	inline DomPoint ConvertToDomain(const RefPoint& refPoint) const
	{
		return this->MeshFace->ConvertToDomain(refPoint);
	}
	inline RefPoint ConvertToReference(const DomPoint& domainPoint) const
	{
		return this->MeshFace->ConvertToReference(domainPoint);
	}

	//-----------//
	//    HHO    //
	//-----------//

public:
	void InitHHO(HHOParameters<Dim>* hho)
	{
		if (this->_massMatrix.rows() > 0)
			return;

		this->HHO = hho;

		//this->ComputeAndSaveQuadraturePoints(basis->GetDegree());
		//this->ComputeAndSaveQuadraturePoints();

		if (hho->OrthonormalizeBases)
		{
			this->Basis = new OrthonormalBasis<Dim - 1>(HHO->FaceBasis, this->MeshFace->Shape());
			//this->_massMatrix = this->MeshFace->Shape()->ComputeMassMatrix(basis);
			//cout << "mass matrix: " << endl << _massMatrix << endl;
			this->_massMatrix = DenseMatrix::Identity(this->Basis->Size(), this->Basis->Size());
			this->_invMassMatrix = DenseMatrix::Identity(this->Basis->Size(), this->Basis->Size());
		}
		else
		{
			this->Basis = hho->FaceBasis;
			this->_massMatrix = this->MeshFace->Shape()->MassMatrix(this->Basis);
			this->_invMassMatrix = this->_massMatrix.inverse();
		}
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
		auto reconstructPhiOnFace = element->EvalPhiOnFace(this->MeshFace, reconstructPhi);

		RefFunction functionToIntegrate = [facePhi, reconstructPhiOnFace](const RefPoint& p) {
			return facePhi->Eval(p) * reconstructPhiOnFace(p);
		};

		int polynomialDegree = facePhi->GetDegree() + reconstructPhi->GetDegree();
		return this->MeshFace->Integral(functionToIntegrate, polynomialDegree);
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

		return this->MeshFace->Integral(functionToIntegrate);
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

	void DeleteUselessMatricesAfterAssembly()
	{
		Utils::Empty(_massMatrix);
	}

	void DeleteUselessMatricesAfterMultigridSetup()
	{
		Utils::Empty(_invMassMatrix);
		this->MeshFace->EmptySavedDomPoints();
	}

private:
	DenseMatrix ComputeTraceMatrix(Element<Dim>* e, FunctionalBasis<Dim>* cellBasis)
	{
		assert(_invMassMatrix.rows() > 0);
		DenseMatrix massCellFace = this->MassMatrix(this->Basis, e, cellBasis);
		return _invMassMatrix * massCellFace;
	}
};