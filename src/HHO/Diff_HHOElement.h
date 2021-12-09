#pragma once
#include "../Mesh/Element.h"
#include "HHOParameters.h"
#include "Diff_HHOFace.h"

template <int Dim>
class Diff_HHOElement
{
public:
	Element<Dim>* MeshElement = nullptr;
	vector<Diff_HHOFace<Dim>*> Faces;

	HHOParameters<Dim>* HHO;
	FunctionalBasis<Dim>* CellBasis;
	FunctionalBasis<Dim>* ReconstructionBasis;

	// Reconstruction operator as a matrix
	DenseMatrix P;

	// Consistency contribution
	DenseMatrix Acons;

	// Stabilization contribution
	DenseMatrix Astab;

	// Local operator matrix (= Acons + Astab)
	DenseMatrix A;

	Eigen::LLT<DenseMatrix> AttSolver;

	Diff_HHOElement() {}

	//-----------------//
	// Element wrapper //
	//-----------------//

	inline BigNumber Number() const
	{
		return this->MeshElement->Number;
	}
	inline int LocalNumberOf(Diff_HHOFace<Dim>* face)
	{
		return this->MeshElement->LocalNumberOf(face->MeshFace);
	}
	inline double SourceTerm(BasisFunction<Dim>* phi, DomFunction f)
	{
		return this->MeshElement->SourceTerm(phi, f);
	}
private:
	inline Tensor<Dim>* DiffTensor() const
	{
		return this->MeshElement->DiffTensor();
	}
	inline DimVector<Dim> OuterNormalVector(Diff_HHOFace<Dim>* face) const
	{
		return this->MeshElement->OuterNormalVector(face->MeshFace);
	}
	inline RefFunction EvalPhiOnFace(Diff_HHOFace<Dim>* face, BasisFunction<Dim>* phi)
	{
		return this->MeshElement->EvalPhiOnFace(face->MeshFace, phi);
	}
	inline function<DimVector<Dim>(RefPoint)> GradPhiOnFace(Diff_HHOFace<Dim>* face, BasisFunction<Dim>* phi)
	{
		return this->MeshElement->GradPhiOnFace(face->MeshFace, phi);
	}

	//----------------------------//
	//   HHO-specific integrals   //
	//----------------------------//

public:
	DenseMatrix MassMatrix(FunctionalBasis<Dim>* basis) const
	{
		return this->MeshElement->MassMatrix(basis);
	}

private:
	double IntegralKGradGradReconstruct(Tensor<Dim>* K, BasisFunction<Dim>* reconstructPhi1, BasisFunction<Dim>* reconstructPhi2) const
	{
		return this->MeshElement->IntegralKGradGradReconstruct(K, reconstructPhi1, reconstructPhi2);
	}

	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis) const
	{
		return this->MeshElement->CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	double ComputeIntegralKGradGrad(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return this->MeshElement->Shape()->ComputeIntegralKGradGrad(K, phi1, phi2);
	}

	//--------------------------------------------------------------------------------//
public:
	void InitHHO(HHOParameters<Dim>* hho, bool assembleLocalMatrix = true)
	{
		this->HHO = hho;

		if (hho->OrthogonalizeElemBases())
		{
			this->ReconstructionBasis = new OrthogonalBasis<Dim>(HHO->ReconstructionBasis, this->MeshElement->Shape(), hho->NElemOrthogonalizations(), hho->OrthonormalizeElemBases());
			this->CellBasis = new FunctionalBasis<Dim>(this->ReconstructionBasis->ExtractLowerBasis(HHO->CellBasis->GetDegree()));
		}
		else
		{
			this->ReconstructionBasis = HHO->ReconstructionBasis;
			this->CellBasis = HHO->CellBasis;
		}
		//DenseMatrix massMatrix = this->MeshElement->Shape()->ComputeMassMatrix(this->ReconstructionBasis);
		//cout << "Reconstruct mass matrix: " << endl << massMatrix << endl;
		//massMatrix = this->MeshElement->Shape()->ComputeMassMatrix(this->CellBasis);
		//cout << "Cell mass matrix: " << endl << massMatrix << endl;


		//this->ComputeAndSaveQuadraturePoints(hho->CellBasis->GetDegree());
		//this->ComputeAndSaveQuadraturePoints(hho->ReconstructionBasis->GetDegree());
		//this->ComputeAndSaveQuadraturePoints();

		if (assembleLocalMatrix)
		{
			this->AssembleReconstructionAndConsistencyMatrices();
			this->AssembleStabilizationMatrix();

			//int nTotalFaceUnknowns = this->Faces.size() * hho->nFaceUnknowns;

			this->A = Acons + Astab;
			auto Att = A.topLeftCorner(hho->nCellUnknowns, hho->nCellUnknowns);
			//auto Aff = A.bottomRightCorner(nTotalFaceUnknowns, nTotalFaceUnknowns);
			//auto Atf = A.topRightCorner(nCellUnknowns, nTotalFaceUnknowns);
			this->AttSolver = Att.llt();
		}
	}

	bool HasOrthogonalBasis() const
	{
		return HHO->OrthogonalizeElemBases() || (this->ReconstructionBasis->IsOrthogonalOnCartesianShapes && this->MeshElement->Shape()->MapsToACartesianShape());
	}

	DenseMatrix SolveCellMassMatrix(const DenseMatrix& M)
	{
		return SolveMassMatrix(this->CellBasis, M);
	}
	DenseMatrix SolveReconstructMassMatrix(const DenseMatrix& M)
	{
		return SolveMassMatrix(this->ReconstructionBasis, M);
	}

private:
	DenseMatrix SolveMassMatrix(FunctionalBasis<Dim>* basis, const DenseMatrix& M)
	{
		if (HHO->OrthonormalizeElemBases())
			return M;
		else if (HHO->OrthogonalizeElemBases())
		{
			Vector d(basis->Size());
			for (BasisFunction<Dim>* phi : basis->LocalFunctions)
				d[phi->LocalNumber] = dynamic_cast<OrthogonalBasisFunction<Dim>*>(phi)->NormSquare;
			return d.asDiagonal().inverse() * M;
		}
		else
			return this->MassMatrix(basis).llt().solve(M);
	}

	Vector ApplyMassMatrix(FunctionalBasis<Dim>* basis, const Vector& v)
	{
		if (HHO->OrthonormalizeElemBases())
			return v;
		else if (HHO->OrthogonalizeElemBases())
		{
			Vector d(basis->Size());
			for (BasisFunction<Dim>* phi : basis->LocalFunctions)
				d[phi->LocalNumber] = dynamic_cast<OrthogonalBasisFunction<Dim>*>(phi)->NormSquare;
			return d.asDiagonal() * v;
		}
		else
			return this->MassMatrix(basis) * v;
	}

public:
	Vector Reconstruct(Vector hybridVector)
	{
		return this->P * hybridVector;
	}
	
	DenseMatrix ReconstructionMatrix()
	{
		return this->P;
	}

	DenseMatrix SolveCellUnknownsMatrix()
	{
		int nTotalFaceUnknowns = this->Faces.size() * HHO->nFaceUnknowns;
		auto Atf = this->A.topRightCorner(HHO->nCellUnknowns, nTotalFaceUnknowns);
		DenseMatrix solveCellUnknowns = -this->AttSolver.solve(Atf);
		return solveCellUnknowns;
	}

	DenseMatrix ReconstructionFromFacesMatrix()
	{
		int nTotalFaceUnknowns = this->Faces.size() * HHO->nFaceUnknowns;

		auto Atf = this->A.topRightCorner(HHO->nCellUnknowns, nTotalFaceUnknowns);
		DenseMatrix solveCellUnknowns = -this->AttSolver.solve(Atf);

		DenseMatrix createHybridVectorFromFacesMatrix(HHO->nCellUnknowns + nTotalFaceUnknowns, nTotalFaceUnknowns);
		createHybridVectorFromFacesMatrix.topRows(HHO->nCellUnknowns) = solveCellUnknowns;
		createHybridVectorFromFacesMatrix.bottomRows(nTotalFaceUnknowns) = DenseMatrix::Identity(nTotalFaceUnknowns, nTotalFaceUnknowns);
		return this->P * createHybridVectorFromFacesMatrix;
	}

	Vector ApplyCellReconstructMassMatrix(const Vector& v)
	{
		if (HHO->OrthonormalizeElemBases())
			return v;
		else if (HHO->OrthogonalizeElemBases())
		{
			Vector d(this->CellBasis->Size());
			for (BasisFunction<Dim>* phi : this->CellBasis->LocalFunctions)
				d[phi->LocalNumber] = dynamic_cast<OrthogonalBasisFunction<Dim>*>(phi)->NormSquare;
			return d.asDiagonal() * v.head(HHO->nCellUnknowns);
		}
		else
		{
			DenseMatrix Nt = this->CellReconstructMassMatrix(this->CellBasis, this->ReconstructionBasis);
			return Nt * v;
		}
	}

	Vector ApplyCellMassMatrix(const Vector& v)
	{
		return this->ApplyMassMatrix(this->CellBasis, v);
	}

private:

	//------------------------------------------------------------------------------------//
	//                Reconstruction operator and consistency contribution                //
	//------------------------------------------------------------------------------------//

	void AssembleReconstructionAndConsistencyMatrices()
	{
		DenseMatrix reconstructionMatrixToInvert = this->AssembleReconstructionMatrixToInvert();
		DenseMatrix rhsMatrix = this->AssembleRHSMatrix();

		DenseMatrix P_and_LagrangeMultipliers = reconstructionMatrixToInvert.colPivHouseholderQr().solve(rhsMatrix);
		// We don't keep the last row, which stores the values of the Lagrange multipliers
		this->P = P_and_LagrangeMultipliers.topRows(HHO->nReconstructUnknowns);

		assert(this->P.cols() == HHO->nCellUnknowns + this->Faces.size() * HHO->nFaceUnknowns);

		DenseMatrix laplacianMatrix = reconstructionMatrixToInvert.topLeftCorner(HHO->nReconstructUnknowns, HHO->nReconstructUnknowns); // Grad-Grad part (block S)

		this->Acons.noalias() = this->P.transpose() * laplacianMatrix * this->P;
	}

	DenseMatrix AssembleReconstructionMatrixToInvert()
	{
		DenseMatrix matrixToInvert(HHO->nReconstructUnknowns + 1, HHO->nReconstructUnknowns + 1);
		this->AssembleGradientReconstructionMatrix(matrixToInvert); // Block S
		this->AssembleMeanValueCondition(matrixToInvert); // Blocks L and L_transpose
		matrixToInvert.bottomRightCorner<1, 1>() << 0;
		return matrixToInvert;
	}

	void AssembleGradientReconstructionMatrix(DenseMatrix & reconstructionMatrixToInvert)
	{
		for (BasisFunction<Dim>* phi1 : this->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : this->ReconstructionBasis->LocalFunctions)
			{
				if (phi2->LocalNumber > phi1->LocalNumber)
					break;
				double value = this->IntegralKGradGradReconstruct(this->DiffTensor(), phi1, phi2);
				reconstructionMatrixToInvert(phi1->LocalNumber, phi2->LocalNumber) = value;
				reconstructionMatrixToInvert(phi2->LocalNumber, phi1->LocalNumber) = value;
			}
		}
	}

	void AssembleMeanValueCondition(DenseMatrix & reconstructionMatrixToInvert)
	{
		int last = this->ReconstructionBasis->Size();
		for (BasisFunction<Dim>* phi : this->ReconstructionBasis->LocalFunctions)
		{
			double meanValue = this->MeshElement->Integral(phi); // TODO: this can be computed only once on the reference element
			reconstructionMatrixToInvert(last, phi->LocalNumber) = meanValue;
			reconstructionMatrixToInvert(phi->LocalNumber, last) = meanValue;
		}
	}

	void AssembleMeanValueConditionRHS(DenseMatrix & rhsMatrix)
	{
		int last = this->ReconstructionBasis->Size();
		for (BasisFunction<Dim>* phi : this->CellBasis->LocalFunctions)
			rhsMatrix(last, phi->LocalNumber) = this->MeshElement->Integral(phi); // TODO: this can be computed only once on the reference element
	}

	DenseMatrix AssembleRHSMatrix()
	{
		auto nTotalFaceUnknowns = this->Faces.size() * HHO->nFaceUnknowns;
		auto nColumns = HHO->nCellUnknowns + nTotalFaceUnknowns;
		DenseMatrix rhsMatrix(HHO->nReconstructUnknowns + 1, nColumns);

		// Top-left corner (Block Bt)
		this->AssembleIntegrationByPartsRHS_cell(rhsMatrix);

		// Top-right corner (Block B_frontier)
		for (auto face : this->Faces)
			this->AssembleIntegrationByPartsRHS_face(rhsMatrix, face);

		// Bottom-left corner (Block Lt)
		this->AssembleMeanValueConditionRHS(rhsMatrix);

		// Bottom-right corner (0)
		rhsMatrix.bottomRightCorner(1, nTotalFaceUnknowns) << Eigen::ArrayXXd::Zero(1, nTotalFaceUnknowns);

		return rhsMatrix;
	}

	void AssembleIntegrationByPartsRHS_cell(DenseMatrix & rhsMatrix)
	{
		for (BasisFunction<Dim>* reconstructPhi : this->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* cellPhi : this->CellBasis->LocalFunctions)
				rhsMatrix(reconstructPhi->LocalNumber, DOFNumber(cellPhi)) = this->IntegrationByPartsRHS_cell(reconstructPhi, cellPhi);
		}
	}

	void AssembleIntegrationByPartsRHS_face(DenseMatrix & rhsMatrix, Diff_HHOFace<Dim> * face)
	{
		for (BasisFunction<Dim>* reconstructPhi : this->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim - 1> * facePhi : face->Basis->LocalFunctions)
				rhsMatrix(reconstructPhi->LocalNumber, DOFNumber(face, facePhi)) = this->IntegrationByPartsRHS_face(face, reconstructPhi, facePhi);
		}
	}

	double IntegrationByPartsRHS_cell(BasisFunction<Dim>* reconstructPhi, BasisFunction<Dim>* cellPhi)
	{
		if (reconstructPhi->GetDegree() == 0)
			return 0;

		double integralGradGrad = this->ComputeIntegralKGradGrad(this->DiffTensor(), reconstructPhi, cellPhi);

		double sumFaces = 0;
		for (auto face : this->Faces)
		{
			auto phi = this->EvalPhiOnFace(face, cellPhi);
			auto gradPhi = this->GradPhiOnFace(face, reconstructPhi);
			auto normal = this->OuterNormalVector(face);

			RefFunction functionToIntegrate = [this, phi, gradPhi, normal](const RefPoint& p) {
				return (this->DiffTensor() * gradPhi(p)).dot(normal) * phi(p);
			};

			int polynomialDegree = reconstructPhi->GetDegree() - 1 + cellPhi->GetDegree();
			double integralFace = face->MeshFace->Integral(functionToIntegrate, polynomialDegree);

			sumFaces += integralFace;
		}

		return integralGradGrad - sumFaces;
	}

	double IntegrationByPartsRHS_face(Diff_HHOFace<Dim>* face, BasisFunction<Dim>* reconstructPhi, BasisFunction<Dim-1>* facePhi)
	{
		if (reconstructPhi->GetDegree() == 0)
			return 0;

		auto gradPhi = this->GradPhiOnFace(face, reconstructPhi);
		auto normal = this->OuterNormalVector(face);

		RefFunction functionToIntegrate = [this, facePhi, gradPhi, normal](const RefPoint& p) {
			return (this->DiffTensor() * gradPhi(p)).dot(normal) * facePhi->Eval(p);
		};

		int polynomialDegree = reconstructPhi->GetDegree() - 1 + facePhi->GetDegree();
		return face->MeshFace->Integral(functionToIntegrate, polynomialDegree);
	}

	//---------------------------------------------//
	//                Stabilization                //
	//---------------------------------------------//

	void AssembleStabilizationMatrix()
	{
		int nHybridUnknowns = HHO->nCellUnknowns + this->Faces.size() * HHO->nFaceUnknowns;

		this->Astab = DenseMatrix::Zero(nHybridUnknowns, nHybridUnknowns);

		if (HHO->Stabilization.compare("hho") == 0)
		{
			DenseMatrix Dt;
			if (!this->HasOrthogonalBasis() || !this->ReconstructionBasis->IsHierarchical)
			{
				DenseMatrix cellMassMatrix = this->MassMatrix(this->CellBasis);
				DenseMatrix Nt = this->CellReconstructMassMatrix(this->CellBasis, this->ReconstructionBasis);

				DenseMatrix ProjT = cellMassMatrix.llt().solve(Nt);
				Dt = ProjT * this->P;
			}
			else
				Dt = this->P.topRows(HHO->nCellUnknowns);

			for (int i = 0; i < Dt.rows(); i++)
				Dt(i, i) -= 1;

			for (auto face : this->Faces)
			{
				auto normal = this->OuterNormalVector(face);
				DenseMatrix Mf = face->MassMatrix();
				DenseMatrix ProjF = face->Trace(this->MeshElement, this->ReconstructionBasis);
				DenseMatrix ProjFT = face->Trace(this->MeshElement, this->CellBasis);

				DenseMatrix Df = ProjF * this->P;
				for (int i = 0; i < Df.rows(); i++)
					Df(i, FirstDOFNumber(face) + i) -= 1;

				DenseMatrix DiffTF = Df - ProjFT * Dt;
				double h = face->Diameter();
				this->Astab += DiffTF.transpose() * Mf * DiffTF * (this->DiffTensor() * normal).dot(normal) / h;
			}
		}
		else if (HHO->Stabilization.compare("hdg") == 0)
		{
			for (auto face : this->Faces)
			{
				auto normal = this->OuterNormalVector(face);
				DenseMatrix Mf = face->MassMatrix();
				DenseMatrix ProjFT = face->Trace(this->MeshElement, this->CellBasis);

				DenseMatrix Fpart = DenseMatrix::Zero(HHO->nFaceUnknowns, nHybridUnknowns);
				Fpart.middleCols(FirstDOFNumber(face), HHO->nFaceUnknowns) = DenseMatrix::Identity(HHO->nFaceUnknowns, HHO->nFaceUnknowns);

				DenseMatrix Tpart = DenseMatrix::Zero(HHO->nCellUnknowns, nHybridUnknowns);
				Tpart.middleCols(0, HHO->nCellUnknowns) = DenseMatrix::Identity(HHO->nCellUnknowns, HHO->nCellUnknowns);

				DenseMatrix DiffTF = Fpart - ProjFT * Tpart;
				double h = face->Diameter();
				this->Astab += DiffTF.transpose() * Mf * DiffTF * (this->DiffTensor() * normal).dot(normal) / h;
			}
		}
		else
			assert(false);
	}

	//-----------------------------------------//
	//                DOFNumber                //
	//-----------------------------------------//

private:
	int DOFNumber(BasisFunction<Dim> * cellPhi)
	{
		return cellPhi->LocalNumber;
	}
	int DOFNumber(Diff_HHOFace<Dim> * face, BasisFunction<Dim - 1> * facePhi)
	{
		return FirstDOFNumber(face) + facePhi->LocalNumber;
	}

public:
	int FirstDOFNumber(Diff_HHOFace<Dim> * face)
	{
		return HHO->nCellUnknowns + this->LocalNumberOf(face) * HHO->nFaceUnknowns;
	}

	//--------------------------------------------------------//
	//        Delete matrices when they become useless        //
	//--------------------------------------------------------//
public:
	void DeleteUselessMatricesAfterAssembly()
	{
		Utils::Empty(Acons);
		Utils::Empty(Astab);
	}

	void DeleteUselessMatricesAfterMultigridSetup()
	{
		Utils::Empty(A);
		this->MeshElement->EmptySavedDomPoints();
		/*if (!needToReconstructSolutionLater)
		{
			Utils::Empty(P);
			Utils::Empty(AttSolver);
		}*/
	}

public:
	~Diff_HHOElement()
	{
		if (HHO->OrthogonalizeElemBases())
		{
			delete ReconstructionBasis;
			CellBasis->LocalFunctions.clear();
			delete CellBasis;
		}
	}
};