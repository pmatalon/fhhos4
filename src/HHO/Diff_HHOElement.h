#pragma once
#include "../Mesh/Element.h"
#include "HHOParameters.h"
#include "Diff_HHOFace.h"

template <int Dim>
class Diff_HHOElement : virtual public Element<Dim>
{
private:
	DenseMatrix _projFromReconstruct;
public:
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

	Diff_HHOElement(BigNumber number) : Element<Dim>(number) {}

	//----------------------------//
	//   HHO-specific integrals   //
	//----------------------------//

private:
	double IntegralKGradGradReconstruct(Tensor<Dim>* K, BasisFunction<Dim>* reconstructPhi1, BasisFunction<Dim>* reconstructPhi2)
	{
		return this->Shape()->IntegralKGradGradReconstruct(K, reconstructPhi1, reconstructPhi2);
	}

	DenseMatrix CellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		return this->Shape()->CellMassMatrix(basis);
	}

	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		return this->Shape()->CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	double ComputeIntegralKGradGrad(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return this->Shape()->ComputeIntegralKGradGrad(K, phi1, phi2);
	}

	//--------------------------------------------------------------------------------//
public:
	void InitHHO(HHOParameters<Dim>* hho)
	{
		this->HHO = hho;

		//this->ComputeAndSaveQuadraturePoints(hho->CellBasis->GetDegree());
		//this->ComputeAndSaveQuadraturePoints(hho->ReconstructionBasis->GetDegree());
		this->ComputeAndSaveQuadraturePoints();

		DenseMatrix cellMassMatrix = this->CellMassMatrix(hho->CellBasis);
		DenseMatrix Nt = this->CellReconstructMassMatrix(hho->CellBasis, hho->ReconstructionBasis);
		this->_projFromReconstruct = cellMassMatrix.inverse() * Nt;

		this->AssembleReconstructionAndConsistencyMatrices();
		this->AssembleStabilizationMatrix();

		Utils::Empty(_projFromReconstruct);

		//int nTotalFaceUnknowns = this->Faces.size() * hho->nFaceUnknowns;

		this->A = Acons + Astab;
		auto Att = A.topLeftCorner(hho->nCellUnknowns, hho->nCellUnknowns);
		//auto Aff = A.bottomRightCorner(nTotalFaceUnknowns, nTotalFaceUnknowns);
		//auto Atf = A.topRightCorner(nCellUnknowns, nTotalFaceUnknowns);
		this->invAtt = Att.inverse();
	}
	
	DenseMatrix ComputeCanonicalInjectionMatrixCoarseToFine(FunctionalBasis<Dim>* cellBasis)
	{
		assert(this->FinerElements.size() > 0);
		DenseMatrix J(cellBasis->Size() * this->FinerElements.size(), cellBasis->Size());

		for (auto e : this->FinerElements)
		{
			Diff_HHOElement<Dim>* fineElement = dynamic_cast<Diff_HHOElement<Dim>*>(e);

			DenseMatrix fineCoarseMass(cellBasis->Size(), cellBasis->Size());
			for (BasisFunction<Dim>* finePhi : cellBasis->LocalFunctions)
			{
				for (BasisFunction<Dim>* coarsePhi : cellBasis->LocalFunctions)
				{
					RefFunction functionToIntegrate = [this, fineElement, finePhi, coarsePhi](const RefPoint& fineRefPoint) {
						DomPoint domPoint = fineElement->ConvertToDomain(fineRefPoint);
						RefPoint coarseRefPoint = this->ConvertToReference(domPoint);
						return finePhi->Eval(fineRefPoint)*coarsePhi->Eval(coarseRefPoint);
					};
					
					int polynomialDegree = finePhi->GetDegree() + coarsePhi->GetDegree();
					double integral = fineElement->Integral(functionToIntegrate, polynomialDegree);
					fineCoarseMass(finePhi->LocalNumber, coarsePhi->LocalNumber) = integral;
				}
			}
			
			DenseMatrix fineMass = fineElement->CellMassMatrix(cellBasis);

			J.block(this->LocalNumberOf(fineElement)*cellBasis->Size(), 0, cellBasis->Size(), cellBasis->Size()) = fineMass.inverse() * fineCoarseMass;
		}

		return J;
	}

	DenseMatrix ComputeL2ProjectionMatrixCoarseToFine(FunctionalBasis<Dim>* cellBasis)
	{
		assert(this->OverlappingFineElements.size() > 0);
		DenseMatrix L2Proj(cellBasis->Size() * this->OverlappingFineElements.size(), cellBasis->Size());

		for (auto e : this->OverlappingFineElements)
		{
			if (e->PhysicalPart != this->PhysicalPart)
			{
				Utils::Warning("This coarse element is overlapped by a fine one that is not in the same physical part. The coarsening/refinement strategy must prevent that.");
				continue;
			}
			Diff_HHOElement<Dim>* fineElement = dynamic_cast<Diff_HHOElement<Dim>*>(e);

			DenseMatrix fineCoarseMass(cellBasis->Size(), cellBasis->Size());
			for (BasisFunction<Dim>* finePhi : cellBasis->LocalFunctions)
			{
				for (BasisFunction<Dim>* coarsePhi : cellBasis->LocalFunctions)
				{
					RefFunction finePhiCoarsePhi = [this, fineElement, finePhi, coarsePhi](const RefPoint& fineRefPoint) {
						DomPoint fineDomPoint = fineElement->ConvertToDomainAndSaveResult(fineRefPoint, true);

						/*if ((fineElement->CoarserElement == this && fineElement->IsFullyEmbeddedInCoarseElement) && !this->Contains(fineDomPoint))
						{
							cout << endl << "% Coarse element: " << endl;
							this->ExportToMatlab("r");
							cout << endl << "% Fine element: " << endl;
							fineElement->ExportToMatlab("b");
							cout << endl << "% Quadrature point: " << endl;
							MatlabScript s;
							s.PlotPoint(fineDomPoint, "k+");
							this->Contains(fineDomPoint);
							this->FullyEmbeds(fineElement);
						}*/
						if ((fineElement->CoarserElement == this && fineElement->IsFullyEmbeddedInCoarseElement) || this->Contains(fineDomPoint))
						{
							RefPoint coarseRefPoint = this->ConvertToReference(fineDomPoint);
							return finePhi->Eval(fineRefPoint)*coarsePhi->Eval(coarseRefPoint);
						}
						else
							return 0.0;
					};

					double integral = 0;
					if (fineElement->CoarserElement == this && fineElement->IsFullyEmbeddedInCoarseElement)
					{
						int degree = finePhi->GetDegree() + coarsePhi->GetDegree();
						integral = fineElement->Integral(finePhiCoarsePhi, degree);
					}
					else
						integral = fineElement->Integral(finePhiCoarsePhi); // use all quadrature points because it's not a polynomial
					fineCoarseMass(finePhi->LocalNumber, coarsePhi->LocalNumber) = integral;
				}
			}

			DenseMatrix fineMass = fineElement->CellMassMatrix(cellBasis);

			L2Proj.block(this->LocalNumberOfOverlapping(fineElement)*cellBasis->Size(), 0, cellBasis->Size(), cellBasis->Size()) = fineMass.inverse() * fineCoarseMass;
		}

		return L2Proj;
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

	inline double MatrixTerm(BasisFunction<Dim>* cellPhi1, BasisFunction<Dim>* cellPhi2)
	{
		return this->A(DOFNumber(cellPhi1), DOFNumber(cellPhi2));
	}
	inline double MatrixTerm(Face<Dim>* face, BasisFunction<Dim>* cellPhi, BasisFunction<Dim - 1> * facePhi)
	{
		return this->A(DOFNumber(cellPhi), DOFNumber(face, facePhi));
	}
	inline double MatrixTerm(Face<Dim>* face1, BasisFunction<Dim - 1> * facePhi1, Face<Dim>* face2, BasisFunction<Dim - 1> * facePhi2)
	{
		return this->A(DOFNumber(face1, facePhi1), DOFNumber(face2, facePhi2));
	}

	inline double ConsistencyTerm(BasisFunction<Dim>* cellPhi1, BasisFunction<Dim>* cellPhi2)
	{
		return this->Acons(DOFNumber(cellPhi1), DOFNumber(cellPhi2));
	}
	inline double ConsistencyTerm(Face<Dim>* face, BasisFunction<Dim>* cellPhi, BasisFunction<Dim - 1> * facePhi)
	{
		return this->Acons(DOFNumber(cellPhi), DOFNumber(face, facePhi));
	}
	inline double ConsistencyTerm(Face<Dim>* face1, BasisFunction<Dim - 1> * facePhi1, Face<Dim>* face2, BasisFunction<Dim - 1> * facePhi2)
	{
		return this->Acons(DOFNumber(face1, facePhi1), DOFNumber(face2, facePhi2));
	}

	inline double StabilizationTerm(BasisFunction<Dim>* cellPhi1, BasisFunction<Dim>* cellPhi2)
	{
		return this->Astab(DOFNumber(cellPhi1), DOFNumber(cellPhi2));
	}
	inline double StabilizationTerm(Face<Dim>* face, BasisFunction<Dim>* cellPhi, BasisFunction<Dim - 1> * facePhi)
	{
		return this->Astab(DOFNumber(cellPhi), DOFNumber(face, facePhi));
	}
	inline double StabilizationTerm(Face<Dim>* face1, BasisFunction<Dim - 1> * facePhi1, Face<Dim>* face2, BasisFunction<Dim - 1> * facePhi2)
	{
		return this->Astab(DOFNumber(face1, facePhi1), DOFNumber(face2, facePhi2));
	}

	inline double ReconstructionTerm(BasisFunction<Dim>* reconstructPhi, BasisFunction<Dim>* cellPhi2)
	{
		return this->P(reconstructPhi->LocalNumber, DOFNumber(cellPhi2));
	}
	inline double ReconstructionTerm(BasisFunction<Dim>* reconstructPhi, Face<Dim>* face, BasisFunction<Dim - 1> * facePhi)
	{
		return this->P(reconstructPhi->LocalNumber, DOFNumber(face, facePhi));
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
				reconstructionMatrixToInvert(phi1->LocalNumber, phi2->LocalNumber) = this->IntegralKGradGradReconstruct(this->DiffTensor(), phi1, phi2);
		}
	}

	void AssembleMeanValueCondition(DenseMatrix & reconstructionMatrixToInvert)
	{
		int last = HHO->ReconstructionBasis->NumberOfLocalFunctionsInElement(NULL);
		for (BasisFunction<Dim>* phi : HHO->ReconstructionBasis->LocalFunctions)
		{
			double meanValue = this->Integral(phi);
			reconstructionMatrixToInvert(last, phi->LocalNumber) = meanValue;
			reconstructionMatrixToInvert(phi->LocalNumber, last) = meanValue;
		}
	}

	void AssembleMeanValueConditionRHS(DenseMatrix & rhsMatrix)
	{
		int last = HHO->ReconstructionBasis->NumberOfLocalFunctionsInElement(NULL);
		for (BasisFunction<Dim>* phi : HHO->CellBasis->LocalFunctions)
			rhsMatrix(last, phi->LocalNumber) = this->Integral(phi);
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

	void AssembleIntegrationByPartsRHS_face(DenseMatrix & rhsMatrix, Face<Dim> * f)
	{
		Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);
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
		for (auto f : this->Faces)
		{
			Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);

			auto phi = this->EvalPhiOnFace(face, cellPhi);
			auto gradPhi = this->GradPhiOnFace(face, reconstructPhi);
			auto normal = this->OuterNormalVector(face);

			RefFunction functionToIntegrate = [this, phi, gradPhi, normal](const RefPoint& p) {
				return (this->DiffTensor() * gradPhi(p)).dot(normal) * phi(p);
			};

			int polynomialDegree = reconstructPhi->GetDegree() - 1 + cellPhi->GetDegree();
			double integralFace = face->Integral(functionToIntegrate, polynomialDegree);

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
		return face->Integral(functionToIntegrate, polynomialDegree);
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
			DenseMatrix ProjT = _projFromReconstruct;
			DenseMatrix Dt = ProjT * this->P;
			for (int i = 0; i < Dt.rows(); i++)
				Dt(i, i) -= 1;

			for (auto f : this->Faces)
			{
				Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);
				auto normal = this->OuterNormalVector(face);
				DenseMatrix Mf = face->FaceMassMatrix();
				DenseMatrix ProjF = face->TraceUsingReconstructBasis(this);
				DenseMatrix ProjFT = face->TraceUsingCellBasis(this);

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
			for (auto f : this->Faces)
			{
				Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);
				auto normal = this->OuterNormalVector(face);
				DenseMatrix Mf = face->FaceMassMatrix();
				DenseMatrix ProjFT = face->TraceUsingCellBasis(this);

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

	//-------------------------------------------------------------------------//
	//  Find polynomial faces which reconstruct the k+1 polynomial on the cell //
	//  and minimize the L2-norm                                               //
	//-------------------------------------------------------------------------//

public:
	DenseMatrix FindFacesPolyWhichReconstructOnTheCell()
	{
		// Assembly of the Lagrangian matrix
		int nBoundaryUnknowns = this->Faces.size() * HHO->nFaceUnknowns;
		DenseMatrix boundaryMassMatrix = DenseMatrix::Zero(nBoundaryUnknowns, nBoundaryUnknowns);
		for (Face<Dim>* face : this->Faces)
		{
			Diff_HHOFace<Dim>* f = dynamic_cast<Diff_HHOFace<Dim>*>(face);
			int localNumber = this->LocalNumberOf(f);
			boundaryMassMatrix.block(localNumber, localNumber, HHO->nFaceUnknowns, HHO->nFaceUnknowns) = f->FaceMassMatrix();
		}

		DenseMatrix C = ReconstructionFromFacesMatrix();

		DenseMatrix lagrangianMatrix(nBoundaryUnknowns + C.rows(), nBoundaryUnknowns + C.rows());
		lagrangianMatrix.topLeftCorner(nBoundaryUnknowns, nBoundaryUnknowns) = boundaryMassMatrix;
		lagrangianMatrix.topRightCorner(C.cols(), C.rows()) = C.transpose();
		lagrangianMatrix.bottomLeftCorner(C.rows(), C.cols()) = C;
		lagrangianMatrix.bottomRightCorner(C.rows(), C.rows()) = DenseMatrix::Zero(C.rows(), C.rows());

		// Assembly of the right-hand side
		DenseMatrix rhs(HHO->nReconstructUnknowns + C.cols(), HHO->nReconstructUnknowns);
		rhs.topRows(C.cols()) = DenseMatrix::Zero(C.cols(), HHO->nReconstructUnknowns);
		rhs.bottomRows(HHO->nReconstructUnknowns) = DenseMatrix::Identity(HHO->nReconstructUnknowns, HHO->nReconstructUnknowns);

		// Solving
		Eigen::ColPivHouseholderQR<DenseMatrix> solver = lagrangianMatrix.colPivHouseholderQr();
		DenseMatrix solutionWithLagrangeCoeffs = solver.solve(rhs);
		DenseMatrix solution = solutionWithLagrangeCoeffs.topRows(nBoundaryUnknowns);

		return solution;
	}

	//-------------------------------------------//
	//  Used by the algorithm from Wildey et al. //
	//-------------------------------------------//

public:
	DenseMatrix StaticallyCondenseInteriorFinerFaces(const SparseMatrix& fineA)
	{
		int nFaceUnknowns = HHO->nFaceUnknowns;

		int nFineInteriorFaces = this->FinerFacesRemoved.size();
		int nFineBoundaryFaces = 0;
		for (Face<Dim>* cf : this->Faces)
			nFineBoundaryFaces += cf->FinerFaces.size();

		DenseMatrix Aii = DenseMatrix::Zero(nFineInteriorFaces*nFaceUnknowns, nFineInteriorFaces*nFaceUnknowns);
		DenseMatrix Aib = DenseMatrix::Zero(nFineInteriorFaces*nFaceUnknowns, nFineBoundaryFaces*nFaceUnknowns);

		for (int i = 0; i < this->FinerFacesRemoved.size(); i++)
		{
			// Aii
			Diff_HHOFace<Dim>* fi = dynamic_cast<Diff_HHOFace<Dim>*>(this->FinerFacesRemoved[i]);
			for (int j = 0; j < this->FinerFacesRemoved.size(); j++)
			{
				Diff_HHOFace<Dim>* fj = dynamic_cast<Diff_HHOFace<Dim>*>(this->FinerFacesRemoved[j]);
				Aii.block(i*nFaceUnknowns, j*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns) = fineA.block(fi->Number*nFaceUnknowns, fj->Number*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns);
			}
			// Aib
			BigNumber fineFaceLocalNumberInCoarseElem = 0;
			for (Face<Dim>* cf : this->Faces)
			{
				if (cf->IsDomainBoundary)
					continue;

				for (int j = 0; j < cf->FinerFaces.size(); j++)
				{
					Diff_HHOFace<Dim>* fj = dynamic_cast<Diff_HHOFace<Dim>*>(cf->FinerFaces[j]);
					Aib.block(i*nFaceUnknowns, fineFaceLocalNumberInCoarseElem*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns) = fineA.block(fi->Number*nFaceUnknowns, fj->Number*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns);
					fineFaceLocalNumberInCoarseElem++;
				}
			}
		}


		int nCoarseFaces = 0;
		for (Face<Dim>* cf : this->Faces)
			nCoarseFaces++;

		DenseMatrix J = DenseMatrix::Zero(nFineBoundaryFaces*nFaceUnknowns, nCoarseFaces*nFaceUnknowns);
		BigNumber fineFaceLocalNumberInCoarseElem = 0;
		for (auto cf : this->Faces)
		{
			if (cf->IsDomainBoundary)
				continue;

			Diff_HHOFace<Dim>* coarseFace = dynamic_cast<Diff_HHOFace<Dim>*>(cf);

			DenseMatrix local_J_f_c = coarseFace->ComputeCanonicalInjectionMatrixCoarseToFine(HHO->FaceBasis);
			for (auto fineFace : coarseFace->FinerFaces)
			{
				BigNumber fineFaceLocalNumberInCoarseFace = coarseFace->LocalNumberOf(fineFace);
				BigNumber coarseFaceLocalNumberInCoarseElem = this->LocalNumberOf(coarseFace);

				J.block(fineFaceLocalNumberInCoarseElem*nFaceUnknowns, coarseFaceLocalNumberInCoarseElem*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns) = local_J_f_c.block(fineFaceLocalNumberInCoarseFace*nFaceUnknowns, 0, nFaceUnknowns, nFaceUnknowns);
				fineFaceLocalNumberInCoarseElem++;
			}
		}

		DenseMatrix resolveCondensedFinerFacesFromFineBoundary = -Aii.inverse()*Aib;
		DenseMatrix resolveCondensedFinerFacesFromCoarseBoundary = resolveCondensedFinerFacesFromFineBoundary * J;

		return resolveCondensedFinerFacesFromCoarseBoundary;
	}

	//-----------------------------------------//
	//                DOFNumber                //
	//-----------------------------------------//

private:
	int DOFNumber(BasisFunction<Dim> * cellPhi)
	{
		return cellPhi->LocalNumber;
	}
	int DOFNumber(Face<Dim> * face, BasisFunction<Dim - 1> * facePhi)
	{
		return FirstDOFNumber(face) + facePhi->LocalNumber;
	}

public:
	int FirstDOFNumber(Face<Dim> * face)
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
		this->EmptySavedDomPoints();
		/*if (!needToReconstructSolutionLater)
		{
			Utils::Empty(P);
			Utils::Empty(invAtt);
		}*/
	}
};