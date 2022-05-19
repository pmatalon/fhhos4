#pragma once
#include "../../Mesh/Element.h"
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
private:
	inline const Tensor<Dim>& DiffTensor() const
	{
		return this->MeshElement->DiffTensor();
	}
	inline DimVector<Dim> OuterNormalVector(Diff_HHOFace<Dim>* face) const
	{
		return this->MeshElement->OuterNormalVector(face->MeshFace);
	}
	inline double EvalTrace(Diff_HHOFace<Dim>* face, BasisFunction<Dim>* phi, const RefPoint& refPointOnFace)
	{
		return this->MeshElement->EvalTrace(face->MeshFace, phi, refPointOnFace);
	}
	inline DimVector<Dim> EvalGradOnFace(Diff_HHOFace<Dim>* face, BasisFunction<Dim>* phi, const RefPoint& refPointOnFace)
	{
		return this->MeshElement->EvalGradOnFace(face->MeshFace, phi, refPointOnFace);
	}
	inline DomPoint ConvertToDomain(const RefPoint& refPoint) const
	{
		return this->MeshElement->ConvertToDomain(refPoint);
	}
	inline RefPoint ConvertToReference(const DomPoint& domainPoint) const
	{
		return this->MeshElement->ConvertToReference(domainPoint);
	}

	//----------------------------//
	//   HHO-specific integrals   //
	//----------------------------//

public:
	DenseMatrix MassMatrix(FunctionalBasis<Dim>* basis) const
	{
		if (HHO->OrthonormalizeElemBases())
			return DenseMatrix::Identity(basis->Size(), basis->Size());
		else if (HHO->OrthogonalizeElemBases())
		{
			Vector d(basis->Size());
			for (BasisFunction<Dim>* phi : basis->LocalFunctions())
				d[phi->LocalNumber] = dynamic_cast<OrthogonalBasisFunction<Dim>*>(phi)->NormSquare;
			return d.asDiagonal();
		}
		else
			return this->MeshElement->MassMatrix(basis);
	}

private:
	double IntegralKGradGrad(const Tensor<Dim>& K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return this->MeshElement->IntegralKGradGrad(K, phi1, phi2);
	}

	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis) const
	{
		return this->MeshElement->CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	//--------------------------------------------------------------------------------//
public:
	void InitHHO(HHOParameters<Dim>* hho, bool assembleLocalMatrix = true)
	{
		this->HHO = hho;

		if (hho->OrthogonalizeElemBases())
		{
			PhysicalShapeWithConstantJacobian<Dim>* shapeCstJac = dynamic_cast<PhysicalShapeWithConstantJacobian<Dim>*>(this->MeshElement->Shape());
			const OrthogonalBasis<Dim>* refShapeOrthogBasis = nullptr;
			if (shapeCstJac)
				refShapeOrthogBasis = shapeCstJac->RefShape()->OrthogonalizedBasis(HHO->ReconstructionBasis);

			/*if (shapeCstJac && refShapeOrthogBasis)
				this->ReconstructionBasis = new OrthogonalBasis<Dim>(*refShapeOrthogBasis, shapeCstJac->DetJacobian(), hho->OrthonormalizeElemBases());
			else*/
				this->ReconstructionBasis = new OrthogonalBasis<Dim>(HHO->ReconstructionBasis, this->MeshElement->Shape(), hho->NElemOrthogonalizations(), hho->OrthonormalizeElemBases());
			this->CellBasis = this->ReconstructionBasis->CreateLowerDegreeBasis(HHO->CellBasis->GetDegree());
		}
		else
		{
			this->ReconstructionBasis = HHO->ReconstructionBasis;
			this->CellBasis = HHO->CellBasis;
		}

		/*cout << "Reconstruct mass matrix (computed): " << endl << this->MeshElement->Shape()->ComputeAndReturnMassMatrix(this->ReconstructionBasis) << endl;
		cout << "Reconstruct mass matrix (applied): " << endl << this->MassMatrix(this->ReconstructionBasis) << endl;
		cout << "Cell mass matrix (computed): " << endl << this->MeshElement->Shape()->ComputeAndReturnMassMatrix(this->CellBasis) << endl;
		cout << "Cell mass matrix (applied): " << endl << this->MassMatrix(this->CellBasis) << endl;*/


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
		return HHO->OrthogonalizeElemBases() || (this->ReconstructionBasis->IsOrthogonalOnCartesianShapes() && this->MeshElement->Shape()->MapsToACartesianShape());
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
			for (BasisFunction<Dim>* phi : basis->LocalFunctions())
				d[phi->LocalNumber] = dynamic_cast<OrthogonalBasisFunction<Dim>*>(phi)->NormSquare;
			return d.asDiagonal().inverse() * M;
		}
		else
			return this->MeshElement->MassMatrix(basis).llt().solve(M);
	}

	Vector ApplyMassMatrix(FunctionalBasis<Dim>* basis, const Vector& v)
	{
		if (HHO->OrthonormalizeElemBases())
			return v;
		else if (HHO->OrthogonalizeElemBases())
		{
			Vector d(basis->Size());
			for (BasisFunction<Dim>* phi : basis->LocalFunctions())
				d[phi->LocalNumber] = dynamic_cast<OrthogonalBasisFunction<Dim>*>(phi)->NormSquare;
			return d.asDiagonal() * v;
		}
		else
			return this->MeshElement->MassMatrix(basis) * v;
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

	DenseMatrix CreateHybridVectorFromFacesMatrix()
	{
		int nTotalFaceUnknowns = this->Faces.size() * HHO->nFaceUnknowns;
		DenseMatrix createHybridVectorFromFacesMatrix(HHO->nCellUnknowns + nTotalFaceUnknowns, nTotalFaceUnknowns);
		createHybridVectorFromFacesMatrix.topRows(HHO->nCellUnknowns) = SolveCellUnknownsMatrix();
		createHybridVectorFromFacesMatrix.bottomRows(nTotalFaceUnknowns) = DenseMatrix::Identity(nTotalFaceUnknowns, nTotalFaceUnknowns);
		return createHybridVectorFromFacesMatrix;
	}

	DenseMatrix ReconstructionFromFacesMatrix()
	{
		return this->P * CreateHybridVectorFromFacesMatrix();
	}

	Vector ReconstructFromFaces(const Vector& faceCoeffs)
	{
		assert(faceCoeffs.rows() == this->Faces.size() * HHO->nFaceUnknowns);

		int nTotalFaceUnknowns = faceCoeffs.rows();
		auto Atf = this->A.topRightCorner(HHO->nCellUnknowns, nTotalFaceUnknowns);

		Vector hybrid(HHO->nCellUnknowns + nTotalFaceUnknowns);
		hybrid.head(HHO->nCellUnknowns) = -this->AttSolver.solve(Atf * faceCoeffs);
		hybrid.tail(nTotalFaceUnknowns) = faceCoeffs;
		return this->P * hybrid;
	}

	Vector ApplyCellReconstructMassMatrix(const Vector& v)
	{
		if (HHO->OrthonormalizeElemBases())
			return v.head(HHO->nCellUnknowns);
		else if (HHO->OrthogonalizeElemBases())
		{
			Vector d(this->CellBasis->Size());
			for (BasisFunction<Dim>* phi : this->CellBasis->LocalFunctions())
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

	Vector ApplyReconstructMassMatrix(const Vector& v)
	{
		return this->ApplyMassMatrix(this->ReconstructionBasis, v);
	}

	double IntegralReconstruct(const Vector& reconstructCoeffs)
	{
		double integral = 0;
		for (BasisFunction<Dim>* phi : this->ReconstructionBasis->LocalFunctions())
			integral += reconstructCoeffs[phi->LocalNumber] * this->MeshElement->Integral(phi); // TODO: this can be computed only once on the reference element
		return integral;
	}


public:
	Vector InnerProductWithBasis(FunctionalBasis<Dim>* basis, DomFunction f)
	{
		return MeshElement->Shape()->InnerProductWithBasis(basis, f);
	}
	Vector ProjectOnReconstructBasis(DomFunction f)
	{
		return SolveReconstructMassMatrix(InnerProductWithBasis(ReconstructionBasis, f));
	}

private:

	//------------------------------------------------------------------------------------//
	//                Reconstruction operator and consistency contribution                //
	//------------------------------------------------------------------------------------//

	void AssembleReconstructionAndConsistencyMatrices()
	{
		DenseMatrix reconstructionMatrixToInvert = AssembleReconstructionMatrixToInvert();
		DenseMatrix rhsMatrix = AssembleRHSMatrix(reconstructionMatrixToInvert);

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
		// Block S (stiffness)
		matrixToInvert.topLeftCorner(HHO->nReconstructUnknowns, HHO->nReconstructUnknowns) << MeshElement->IntegralKGradGradMatrix(DiffTensor(), ReconstructionBasis);
		// Blocks L and L_transpose (mean values)
		Vector meanValues = MeshElement->Integral(ReconstructionBasis);
		matrixToInvert.topRightCorner(meanValues.rows(), 1) = meanValues;
		matrixToInvert.bottomLeftCorner(1, meanValues.rows()) = meanValues.transpose();
		// 0
		matrixToInvert.bottomRightCorner<1, 1>() << 0;
		return matrixToInvert;
	}

	DenseMatrix AssembleRHSMatrix(DenseMatrix& reconstructionMatrixToInvert)
	{
		auto nTotalFaceUnknowns = this->Faces.size() * HHO->nFaceUnknowns;
		DenseMatrix rhsMatrix(HHO->nReconstructUnknowns + 1, HHO->nCellUnknowns + nTotalFaceUnknowns);

		// Top-left corner (Block Bt)
		AssembleIntegrationByPartsRHS_cell(rhsMatrix, reconstructionMatrixToInvert);

		// Top-right corner (Block B_frontier)
		for (auto face : this->Faces)
			AssembleIntegrationByPartsRHS_face(rhsMatrix, face);

		// Bottom-left corner (Block Lt) Mean value condition
		if (ReconstructionBasis->IsHierarchical())
			rhsMatrix.bottomLeftCorner(1, CellBasis->Size()) = reconstructionMatrixToInvert.bottomLeftCorner(1, CellBasis->Size());
		else
			rhsMatrix.bottomLeftCorner(1, CellBasis->Size()) = MeshElement->Integral(CellBasis).transpose();

		// Bottom-right corner (0)
		rhsMatrix.bottomRightCorner(1, nTotalFaceUnknowns) << Eigen::ArrayXXd::Zero(1, nTotalFaceUnknowns);

		return rhsMatrix;
	}

	void AssembleIntegrationByPartsRHS_cell(DenseMatrix& rhsMatrix, DenseMatrix& reconstructionMatrixToInvert)
	{
		for (BasisFunction<Dim>* reconstructPhi : this->ReconstructionBasis->LocalFunctions())
		{
			for (BasisFunction<Dim>* cellPhi : this->CellBasis->LocalFunctions())
				rhsMatrix(reconstructPhi->LocalNumber, cellPhi->LocalNumber) = IntegrationByPartsRHS_cell(reconstructPhi, cellPhi, reconstructionMatrixToInvert);
		}
	}

	void AssembleIntegrationByPartsRHS_face(DenseMatrix& rhsMatrix, Diff_HHOFace<Dim>* face)
	{
		for (BasisFunction<Dim>* reconstructPhi : this->ReconstructionBasis->LocalFunctions())
		{
			for (BasisFunction<Dim - 1> * facePhi : face->Basis->LocalFunctions())
				rhsMatrix(reconstructPhi->LocalNumber, DOFNumber(face, facePhi)) = IntegrationByPartsRHS_face(face, reconstructPhi, facePhi);
		}
	}

	double IntegrationByPartsRHS_cell(BasisFunction<Dim>* reconstructPhi, BasisFunction<Dim>* cellPhi, DenseMatrix& reconstructionMatrixToInvert)
	{
		if (reconstructPhi->GetDegree() == 0)
			return 0;

		double integralGradGrad;
		if (ReconstructionBasis->IsHierarchical())
			integralGradGrad = reconstructionMatrixToInvert(reconstructPhi->LocalNumber, cellPhi->LocalNumber);
		else
			integralGradGrad = this->IntegralKGradGrad(this->DiffTensor(), reconstructPhi, cellPhi);

		double sumFaces = 0;
		for (auto face : this->Faces)
		{
			auto normal = this->OuterNormalVector(face);

			RefFunction functionToIntegrate = [this, face, cellPhi, reconstructPhi, &normal](const RefPoint& p) {
				return (this->DiffTensor() * this->EvalGradOnFace(face, reconstructPhi, p)).dot(normal) * this->EvalTrace(face, cellPhi, p);
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

		auto normal = this->OuterNormalVector(face);

		RefFunction functionToIntegrate = [this, face, facePhi, reconstructPhi, &normal](const RefPoint& p) {
			return (this->DiffTensor() * this->EvalGradOnFace(face, reconstructPhi, p)).dot(normal) * facePhi->Eval(p);
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
			if (!this->HasOrthogonalBasis() || !this->ReconstructionBasis->IsHierarchical())
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
				DenseMatrix ProjFT;
				if (ReconstructionBasis->IsHierarchical())
					ProjFT = ProjF.leftCols(CellBasis->Size());
				else
					ProjFT = face->Trace(this->MeshElement, this->CellBasis);

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
			delete CellBasis;
		}
	}
};