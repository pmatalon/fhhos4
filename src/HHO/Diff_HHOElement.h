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

	// Reconstruction operator as a matrix
	DenseMatrix P;

	// Consistency contribution
	DenseMatrix Acons;

	// Stabilization contribution
	DenseMatrix Astab;

	// Local operator matrix (= Acons + Astab)
	DenseMatrix A;

	DenseMatrix invAtt;

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
	void InitHHO(HHOParameters<Dim>* hho)
	{
		this->HHO = hho;

		//this->ComputeAndSaveQuadraturePoints(hho->CellBasis->GetDegree());
		//this->ComputeAndSaveQuadraturePoints(hho->ReconstructionBasis->GetDegree());
		//this->ComputeAndSaveQuadraturePoints();

		this->AssembleReconstructionAndConsistencyMatrices();
		this->AssembleStabilizationMatrix();

		//int nTotalFaceUnknowns = this->Faces.size() * hho->nFaceUnknowns;

		this->A = Acons + Astab;
		auto Att = A.topLeftCorner(hho->nCellUnknowns, hho->nCellUnknowns);
		//auto Aff = A.bottomRightCorner(nTotalFaceUnknowns, nTotalFaceUnknowns);
		//auto Atf = A.topRightCorner(nCellUnknowns, nTotalFaceUnknowns);
		this->invAtt = Att.inverse();
	}
	
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
		DenseMatrix solveCellUnknowns = -this->invAtt * Atf;
		return solveCellUnknowns;
	}

	DenseMatrix ReconstructionFromFacesMatrix()
	{
		int nTotalFaceUnknowns = this->Faces.size() * HHO->nFaceUnknowns;

		auto Atf = this->A.topRightCorner(HHO->nCellUnknowns, nTotalFaceUnknowns);
		DenseMatrix solveCellUnknowns = -this->invAtt * Atf;

		DenseMatrix createHybridVectorFromFacesMatrix(HHO->nCellUnknowns + nTotalFaceUnknowns, nTotalFaceUnknowns);
		createHybridVectorFromFacesMatrix.topRows(HHO->nCellUnknowns) = solveCellUnknowns;
		createHybridVectorFromFacesMatrix.bottomRows(nTotalFaceUnknowns) = DenseMatrix::Identity(nTotalFaceUnknowns, nTotalFaceUnknowns);
		return this->P * createHybridVectorFromFacesMatrix;
	}

	virtual ~Diff_HHOElement() {}

private:

	//------------------------------------------------------------------------------------//
	//                Reconstruction operator and consistency contribution                //
	//------------------------------------------------------------------------------------//

	void AssembleReconstructionAndConsistencyMatrices()
	{
		this->P = DenseMatrix(HHO->nReconstructUnknowns, HHO->nCellUnknowns + this->Faces.size() * HHO->nFaceUnknowns);

		DenseMatrix reconstructionMatrixToInvert = this->AssembleReconstructionMatrixToInvert();
		DenseMatrix rhsMatrix = this->AssembleRHSMatrix();

		Eigen::ColPivHouseholderQR<DenseMatrix> solver = reconstructionMatrixToInvert.colPivHouseholderQr();
		for (int j = 0; j < rhsMatrix.cols(); j++)
		{
			Vector col = solver.solve(rhsMatrix.col(j));

			// We don't keep the last element of the column, which is the Lagrange multiplier
			auto colWithoutLagrangeMultiplier = col.head(HHO->nReconstructUnknowns);
			this->P.col(j) = colWithoutLagrangeMultiplier;
		}

		DenseMatrix laplacianMatrix = reconstructionMatrixToInvert.topLeftCorner(HHO->nReconstructUnknowns, HHO->nReconstructUnknowns); // Grad-Grad part (block S)

		this->Acons = this->P.transpose() * laplacianMatrix * this->P;
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
		for (BasisFunction<Dim>* phi1 : HHO->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : HHO->ReconstructionBasis->LocalFunctions)
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
		int last = HHO->ReconstructionBasis->Size();
		for (BasisFunction<Dim>* phi : HHO->ReconstructionBasis->LocalFunctions)
		{
			double meanValue = this->MeshElement->Integral(phi); // TODO: this can be computed only once on the reference element
			reconstructionMatrixToInvert(last, phi->LocalNumber) = meanValue;
			reconstructionMatrixToInvert(phi->LocalNumber, last) = meanValue;
		}
	}

	void AssembleMeanValueConditionRHS(DenseMatrix & rhsMatrix)
	{
		int last = HHO->ReconstructionBasis->Size();
		for (BasisFunction<Dim>* phi : HHO->CellBasis->LocalFunctions)
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
		for (BasisFunction<Dim>* reconstructPhi : HHO->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* cellPhi : HHO->CellBasis->LocalFunctions)
				rhsMatrix(reconstructPhi->LocalNumber, DOFNumber(cellPhi)) = this->IntegrationByPartsRHS_cell(reconstructPhi, cellPhi);
		}
	}

	void AssembleIntegrationByPartsRHS_face(DenseMatrix & rhsMatrix, Diff_HHOFace<Dim> * face)
	{
		for (BasisFunction<Dim>* reconstructPhi : HHO->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim - 1> * facePhi : HHO->FaceBasis->LocalFunctions)
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
			DenseMatrix cellMassMatrix = this->MassMatrix(HHO->CellBasis);
			DenseMatrix Nt = this->CellReconstructMassMatrix(HHO->CellBasis, HHO->ReconstructionBasis);

			DenseMatrix ProjT = cellMassMatrix.inverse() * Nt;
			DenseMatrix Dt = ProjT * this->P;
			for (int i = 0; i < Dt.rows(); i++)
				Dt(i, i) -= 1;

			for (auto face : this->Faces)
			{
				auto normal = this->OuterNormalVector(face);
				DenseMatrix Mf = face->MassMatrix();
				DenseMatrix ProjF = face->TraceUsingReconstructBasis(this->MeshElement);
				DenseMatrix ProjFT = face->TraceUsingCellBasis(this->MeshElement);

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
				DenseMatrix ProjFT = face->TraceUsingCellBasis(this->MeshElement);

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
			Utils::Empty(invAtt);
		}*/
	}
};