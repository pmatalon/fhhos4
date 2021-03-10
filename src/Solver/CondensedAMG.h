#pragma once
#include <mutex>
#include "HybridAlgebraicMesh.h"
#include "AlgebraicMesh.h"
#include "Multigrid.h"
using namespace std;

class CondensedLevel : public Level
{
private:
	CAMGFaceProlongation _faceProlong = CAMGFaceProlongation::FaceAggregates;
	CAMGProlongation _multigridProlong = CAMGProlongation::ReconstructionTrace1Step;
	int _cellBlockSize;
	int _faceBlockSize;
	double _strongCouplingThreshold;
public:
	const SparseMatrix* A_T_T;
	const SparseMatrix* A_T_F;
	const SparseMatrix* A_F_F;
	const SparseMatrix* inv_A_T_T;
private:
	SparseMatrix Q_T;
	SparseMatrix Q_F;

	SparseMatrix A_T_Tc;
	SparseMatrix A_T_Fc;
	SparseMatrix A_F_Fc;
	SparseMatrix* inv_A_T_Tc;
public:
	SparseMatrix Ac;

public:
	CondensedLevel(int number, int cellBlockSize, int faceBlockSize, double strongCouplingThreshold, CAMGFaceProlongation faceProlong, CAMGProlongation prolongation)
		: Level(number)
	{
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
		this->_multigridProlong = prolongation;
		this->_faceProlong = faceProlong;
	}

	BigNumber NUnknowns() override
	{
		return A_F_F->rows();
	}

	void ExportVector(const Vector& v, string suffix, int levelNumber) override
	{ }

	void ExportMatrix(const SparseMatrix& M, string suffix, int levelNumber) override
	{ 
		Eigen::saveMarket(M, Utils::ProgramArgs.OutputDirectory + "/" + suffix + ".dat");
	}

	void CoarsenMesh(CoarseningStrategy coarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, double requestedCoarseningFactor, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) override
	{
		noCoarserMeshProvided = false;
		coarsestPossibleMeshReached = false;

		if (!this->OperatorMatrix)
		{
			if (this->UseGalerkinOperator)
				ComputeGalerkinOperator();
			else
				SetupDiscretizedOperator();
		}

		const SparseMatrix* A_T_T1 = A_T_T;
		const SparseMatrix* A_T_F1 = A_T_F;
		const SparseMatrix* A_F_F1 = A_F_F;
		SparseMatrix *A_T_T2, *A_T_F2, *A_F_F2, *P, *Q_F, *schur2;
		const SparseMatrix* schur1 = this->OperatorMatrix;
		double actualCoarseningFactor = 0;
		double nCoarsenings = 0;

		while (!CoarseningCriteriaReached(coarseningStgy, requestedCoarseningFactor, nCoarsenings, actualCoarseningFactor))
		{
			std::tie(A_T_T2, A_T_F2, A_F_F2, P, Q_F, schur2, coarsestPossibleMeshReached) = Coarsen(*A_T_T1, *A_T_F1, *A_F_F1, *schur1, coarseningStgy, faceCoarseningStgy);
			if (coarsestPossibleMeshReached)
				return;

			if (nCoarsenings > 0)
				delete A_T_T1, A_T_F1, A_F_F1, schur1;

			A_T_T1 = A_T_T2;
			A_T_F1 = A_T_F2;
			A_F_F1 = A_F_F2;
			schur1 = schur2;

			if (nCoarsenings == 0)
			{
				this->P = *P;
				this->Q_F = *Q_F;
			}
			else
			{
				this->P = this->P * *P;
				this->Q_F = this->Q_F * *Q_F;
			}

			/*if (this->_faceProlong == CAMGFaceProlongation::FaceAggregates || this->_multigridProlong == CAMGProlongation::ReconstructSmoothedTraceOrInject)
			{
				const SparseMatrix* oldSchur = schur;
				schur = new SparseMatrix(P->transpose() * (*oldSchur) * (*P));
				if (nCoarsenings > 0)
					delete oldSchur;
			}*/

			delete P, Q_F;

			if (A_T_F1)
			{
				double nFine = this->A_T_F->cols();
				double nCoarse = A_T_F1->cols();
				actualCoarseningFactor = nFine / nCoarse;
			}
			nCoarsenings++;
		}

		this->A_T_Tc = *A_T_T2;
		this->A_T_Fc = *A_T_F2;
		this->A_F_Fc = *A_F_F2;
		this->Ac = SparseMatrix(*schur2);
		delete schur2;
	}


	bool CoarseningCriteriaReached(CoarseningStrategy coarseningStgy, double requestedCoarseningRatio, int nCoarseningsPerformed, double coarseningRatio)
	{
		if (coarseningStgy == CoarseningStrategy::DoublePairwiseAggregation)
			return nCoarseningsPerformed == 2;
		if (coarseningStgy == CoarseningStrategy::MultiplePairwiseAggregation)
			return coarseningRatio >= requestedCoarseningRatio;
		if (coarseningStgy == CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours)
			return nCoarseningsPerformed == 1; // only 1 pass of coarsening
		return coarseningRatio >= requestedCoarseningRatio;
	}




	// Returns <A_T_Tc, A_T_Fc, A_F_Fc, P, Q_F, schurc, coarsestPossibleMeshReached>
	tuple<SparseMatrix*, SparseMatrix*, SparseMatrix*, SparseMatrix*, SparseMatrix*, SparseMatrix*, bool> Coarsen(const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F, const SparseMatrix& schur, CoarseningStrategy elemCoarseningStgy, FaceCoarseningStrategy faceCoarseningStgy)
	{
		//ExportMatrix(A_T_T, "A_T_T", 0);
		//ExportMatrix(A_T_F, "A_T_F", 0);
		//ExportMatrix(A_F_F, "A_F_F", 0);

		bool coarsestPossibleMeshReached = false;

		bool onlyFacesUsed = this->_faceProlong == CAMGFaceProlongation::FaceAggregates && _multigridProlong == CAMGProlongation::FaceProlongation;

		HybridAlgebraicMesh mesh(_cellBlockSize, _faceBlockSize, _strongCouplingThreshold);
		if (!onlyFacesUsed)
		{
			mesh.Build(A_T_T, A_T_F, A_F_F);

			mesh.Coarsen(elemCoarseningStgy, faceCoarseningStgy, coarsestPossibleMeshReached);
			if (coarsestPossibleMeshReached)
				return { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, coarsestPossibleMeshReached };
		}

		/*cout << "------------- Elem" << endl;
		for (HybridElementAggregate& agg : mesh.CoarseElements)
		{
			cout << "(";
			for (auto e : agg.FineElements)
				cout << e->Number << ", ";
			cout << ")" << endl;
		}

		cout << "------------- Faces" << endl;
		for (HybridFaceAggregate& agg : mesh.CoarseFaces)
		{
			cout << "(";
			for (auto f : agg.FineFaces)
				cout << f->Number << ", ";
			cout << ")" << endl;
		}*/

		// Cell-prolongation operator with only one 1 coefficient per row
		SparseMatrix Q_T = BuildQ_T(mesh);

		SparseMatrix* Q_F;
		if (this->_faceProlong == CAMGFaceProlongation::BoundaryAggregatesInteriorAverage)
		{
			// Face-prolongation operator with only one 1 coefficient per row for kept or aggregated faces, average for removed faces
			Q_F = new SparseMatrix(BuildQ_F(mesh, A_F_F));
		}
		else if (this->_faceProlong == CAMGFaceProlongation::BoundaryAggregatesInteriorZero)
		{
			Q_F = new SparseMatrix(BuildQ_F_0Interior(mesh));
		}
		else if (this->_faceProlong == CAMGFaceProlongation::FaceAggregates)
		{
			AlgebraicMesh skeleton(_faceBlockSize, 0);
			//skeleton.Build(*A_F_F);
			skeleton.Build(schur);
			skeleton.PairWiseAggregate(coarsestPossibleMeshReached);
			if (coarsestPossibleMeshReached)
				return { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, coarsestPossibleMeshReached };
			Q_F = new SparseMatrix(BuildQ_F_AllAggregated(skeleton));
		}
		else
			Utils::FatalError("Unmanaged -face-prolong");

		// Intermediate coarse operators
		SparseMatrix* A_T_Tc = nullptr;
		SparseMatrix* A_T_Fc_tmp = nullptr;
		if (!onlyFacesUsed)
		{
			A_T_Tc     = new SparseMatrix(Q_T.transpose() * A_T_T.selfadjointView<Eigen::Lower>() * Q_T);
			A_T_Fc_tmp = new SparseMatrix(Q_T.transpose() * A_T_F * *Q_F);
		}

		/*ExportMatrix(Q_T, "Q_T", 0);
		ExportMatrix(*Q_F, "Q_F", 0);
		ExportMatrix(*A_T_Tc, "A_T_Tc", 0);
		ExportMatrix(*A_T_Fc_tmp, "A_T_Fc", 0);*/

		SparseMatrix* P;
		if (_multigridProlong == CAMGProlongation::ReconstructionTrace1Step || _multigridProlong == CAMGProlongation::ReconstructionTrace2Steps) // 1
		{
			SparseMatrix inv_A_T_Tc = Utils::InvertBlockDiagMatrix(*A_T_Tc, _cellBlockSize);
			// Theta: reconstruction from the coarse faces to the coarse cells
			SparseMatrix Theta = -inv_A_T_Tc * *A_T_Fc_tmp;
			// Pi: average on both sides of each face
			SparseMatrix Pi = BuildTrace(mesh);

			P = new SparseMatrix(Pi * Q_T * Theta);
		}
		else if (_multigridProlong == CAMGProlongation::FaceProlongation) // 3
		{
			// -g 1 -prolong 3 -face-prolong 3 -cs z
			P = Q_F;
		}
		else if (_multigridProlong == CAMGProlongation::ReconstructionTranspose2Steps) // 4
		{
			/*SparseMatrix inv_A_T_T1 = Utils::InvertBlockDiagMatrix(A_T_T1, _cellBlockSize);
			// Theta: reconstruction from the coarse faces to the coarse cells
			SparseMatrix Theta1 = -inv_A_T_T1 * A_T_F1;
			// Pi: transpose of Theta
			SparseMatrix inv_A_T_T = Utils::InvertBlockDiagMatrix(*A_T_T, _cellBlockSize);
			SparseMatrix Theta = -inv_A_T_T * *A_T_F;
			SparseMatrix P1 = Theta.transpose() * Q_T1 * Theta1;

			this->inv_A_T_Tc = new SparseMatrix(Utils::InvertBlockDiagMatrix(A_T_Tc, _cellBlockSize));
			SparseMatrix Theta2 = -(*inv_A_T_Tc) * A_T_Fc;
			SparseMatrix P2 = Theta1.transpose() * Q_T2 * Theta2;
			this->P = P1 * P2;*/
		}
		else if (_multigridProlong == CAMGProlongation::ReconstructTraceOrInject) // 5
		{
			SparseMatrix inv_A_T_Tc = Utils::InvertBlockDiagMatrix(*A_T_Tc, _cellBlockSize);
			SparseMatrix Theta = -inv_A_T_Tc * *A_T_Fc_tmp;   // Reconstruct
			SparseMatrix Pi = BuildCoarseTraceOnFineRemovedFaces(mesh); // Trace

			SparseMatrix ReconstructAndTrace = Pi * Theta;

			NumberParallelLoop<CoeffsChunk> parallelLoop(mesh.Faces.size());
			parallelLoop.ReserveChunkCoeffsSize(_faceBlockSize * 2 * _faceBlockSize);
			parallelLoop.Execute([this, &ReconstructAndTrace, &Q_F, &mesh](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
				{
					const HybridAlgebraicFace* face = &mesh.Faces[faceNumber];
					if (face->IsRemovedOnCoarseMesh)
						chunk->Results.Coeffs.CopyRows(faceNumber*_faceBlockSize, _faceBlockSize, ReconstructAndTrace);
					else
						chunk->Results.Coeffs.CopyRows(faceNumber*_faceBlockSize, _faceBlockSize, *Q_F);
				});
			P = new SparseMatrix(Q_F->rows(), Q_F->cols());
			parallelLoop.Fill(*P);
		}
		else if (_multigridProlong == CAMGProlongation::ReconstructSmoothedTraceOrInject) // 6
		{
			SparseMatrix inv_A_T_Tc = Utils::InvertBlockDiagMatrix(*A_T_Tc, _cellBlockSize);
			SparseMatrix Theta = -inv_A_T_Tc * *A_T_Fc_tmp;   // Reconstruct
			SparseMatrix Pi = BuildCoarseTraceOnFineRemovedFaces(mesh); // Trace
			//SparseMatrix Pi = BuildCoarseTraceOnFineFaces(mesh); // Trace

			SparseMatrix ReconstructAndTrace = Pi * Theta;

			NumberParallelLoop<CoeffsChunk> parallelLoop(mesh.Faces.size());
			parallelLoop.ReserveChunkCoeffsSize(_faceBlockSize * 2 * _faceBlockSize);
			parallelLoop.Execute([this, &ReconstructAndTrace, &Q_F, &mesh](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
				{
					const HybridAlgebraicFace* face = &mesh.Faces[faceNumber];
					if (face->IsRemovedOnCoarseMesh)
						chunk->Results.Coeffs.CopyRows(faceNumber*_faceBlockSize, _faceBlockSize, ReconstructAndTrace);
					else
						chunk->Results.Coeffs.CopyRows(faceNumber*_faceBlockSize, _faceBlockSize, *Q_F);
				});
			SparseMatrix ReconstructTraceOrInject(Q_F->rows(), Q_F->cols());
			parallelLoop.Fill(ReconstructTraceOrInject);

			BlockJacobi blockJacobi(_faceBlockSize, 2.0/3.0);
			blockJacobi.Setup(schur);
			SparseMatrix J = blockJacobi.IterationMatrix(); // Smoothing

			SparseMatrix ReconstructAndSmoothedTrace = J * ReconstructTraceOrInject;

			NumberParallelLoop<CoeffsChunk> parallelLoop2(mesh.Faces.size());
			parallelLoop2.ReserveChunkCoeffsSize(_faceBlockSize * 2 * _faceBlockSize);
			parallelLoop2.Execute([this, &ReconstructAndSmoothedTrace, &Q_F, &mesh](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
				{
					const HybridAlgebraicFace* face = &mesh.Faces[faceNumber];
					if (face->IsRemovedOnCoarseMesh)
						chunk->Results.Coeffs.CopyRows(faceNumber*_faceBlockSize, _faceBlockSize, ReconstructAndSmoothedTrace);
					else
						chunk->Results.Coeffs.CopyRows(faceNumber*_faceBlockSize, _faceBlockSize, *Q_F);
				});
			P = new SparseMatrix(Q_F->rows(), Q_F->cols());
			parallelLoop2.Fill(*P);
		}
		else if (_multigridProlong == CAMGProlongation::FindInteriorThatReconstructs) // 7
		{
			int cbs = _cellBlockSize;
			int fbs = _faceBlockSize;

			SparseMatrix inv_A_T_Tc = Utils::InvertBlockDiagMatrix(*A_T_Tc, cbs);
			SparseMatrix inv_A_T_T  = Utils::InvertBlockDiagMatrix(A_T_T, cbs);

			NumberParallelLoop<CoeffsChunk> parallelLoop(mesh.CoarseElements.size());
			//parallelLoop.ReserveChunkCoeffsSize(fbs * 2 * fbs);
			parallelLoop.Execute([this, cbs, fbs, &A_T_F, &Q_F, &Q_T, &inv_A_T_T, &inv_A_T_Tc, A_T_Fc_tmp, &mesh](BigNumber ceNumber, ParallelChunk<CoeffsChunk>* chunk)
				{
					const HybridElementAggregate& ce = mesh.CoarseElements[ceNumber];

					if (ce.RemovedFineFaces.empty())
						return;

					// Construction of Theta_Tc
					DenseMatrix invA_Tc_Tc = inv_A_T_Tc.block(ce.Number, ce.Number, cbs, cbs);

					DenseMatrix A_Tc_F(cbs, ce.CoarseFaces.size()*fbs);
					for (int cfLocalNumber = 0; cfLocalNumber < ce.CoarseFaces.size(); cfLocalNumber++)
					{
						HybridFaceAggregate* cf = ce.CoarseFaces[cfLocalNumber];
						A_Tc_F.block(0, cfLocalNumber*fbs, cbs, fbs) = A_T_Fc_tmp->block(ce.Number*cbs, cf->Number*fbs, cbs, fbs);
					}
					DenseMatrix Theta_Tc = -invA_Tc_Tc * A_Tc_F;

					// Matrix for minimization problem
					DenseMatrix M(ce.FineElements.size()*cbs, ce.RemovedFineFaces.size()*fbs);

					// RHS for minimization problem
					DenseMatrix f = DenseMatrix::Zero(ce.FineElements.size()*cbs, ce.CoarseFaces.size()*fbs);

					for (int feLocalNumber = 0; feLocalNumber < ce.FineElements.size(); ++feLocalNumber)
					{
						HybridAlgebraicElement* fe = ce.FineElements[feLocalNumber];

						// Construction of Theta_Tf
						DenseMatrix invA_Tf_Tf = inv_A_T_T.block(fe->Number, fe->Number, cbs, cbs);
						DenseMatrix A_Tf_F(cbs, fe->Faces.size()*fbs);
						for (int ffLocalNumber = 0; ffLocalNumber < fe->Faces.size(); ffLocalNumber++)
						{
							HybridAlgebraicFace* ff = fe->Faces[ffLocalNumber];
							A_Tf_F.block(0, ffLocalNumber*fbs, cbs, fbs) = A_T_F.block(fe->Number*cbs, ff->Number*fbs, cbs, fbs);
						}
						DenseMatrix Theta_Tf = -invA_Tf_Tf * A_Tf_F;
						//cout << "Theta_Tf = " << endl << Theta_Tf << endl;
						assert((Theta_Tf.array() > 0).all());

						// Construction of Theta_Tf_int (part of Theta_Tf of interior faces)
						//             and Theta_Tf_ext (part of Theta_Tf of exterior faces)
						int nInteriorFaces = 0;
						int nExteriorFaces = 0;
						for (HybridAlgebraicFace* ff : fe->Faces)
						{
							if (ff->IsRemovedOnCoarseMesh)
								nInteriorFaces++;
							else
								nExteriorFaces++;
						}
						DenseMatrix Theta_Tf_int = DenseMatrix::Zero(cbs, ce.RemovedFineFaces.size()*fbs);
						DenseMatrix Theta_Tf_ext = DenseMatrix::Zero(cbs, nExteriorFaces*fbs);
						int localExtFFNumber = 0;
						DenseMatrix Q_F_restrict_partialTc_corestrict_partialTf = DenseMatrix::Zero(nExteriorFaces*fbs, ce.CoarseFaces.size()*fbs);
						for (int ffLocalNumber = 0; ffLocalNumber < fe->Faces.size(); ffLocalNumber++)
						{
							HybridAlgebraicFace* ff = fe->Faces[ffLocalNumber];
							if (ff->IsRemovedOnCoarseMesh)
							{
								int localNumberInCE = ce.LocalRemovedFineFaceNumber(ff);
								Theta_Tf_int.block(0, localNumberInCE*fbs, cbs, fbs) = Theta_Tf.block(0, ffLocalNumber*fbs, cbs, fbs);
								assert(Theta_Tf_int.norm() != 0);
							}
							else
							{
								//int localNumberInCE = ce.LocalFineFaceNumber(ff);
								Theta_Tf_ext.block(0, localExtFFNumber*fbs, cbs, fbs) = Theta_Tf.block(0, ffLocalNumber*fbs, cbs, fbs);

								Q_F_restrict_partialTc_corestrict_partialTf.block(localExtFFNumber*fbs, ce.LocalCoarseFaceNumber(ff->CoarseFace)*fbs, fbs, fbs) = Q_F->block(ff->Number*fbs, ff->CoarseFace->Number*fbs, fbs, fbs);
								localExtFFNumber++;
							}
						}

						// Part of matrix for minimization problem
						M.middleRows(feLocalNumber*cbs, cbs) = Theta_Tf_int;

						// Part of RHS for minimization problem
						DenseMatrix Q_Tc_corestrict_Tf = Q_T.block(fe->Number*cbs, ce.Number*cbs, cbs, cbs);
						DenseMatrix coarseReconstructionThenInjection = Q_Tc_corestrict_Tf * Theta_Tc;
						DenseMatrix faceProlongThenFineReconstruction = Theta_Tf_ext * Q_F_restrict_partialTc_corestrict_partialTf;
						f.middleRows(feLocalNumber*cbs, cbs) = coarseReconstructionThenInjection - faceProlongThenFineReconstruction;
					}

					DenseMatrix x = M.colPivHouseholderQr().solve(f);
					/*if (std::isnan(x.norm()) || std::isinf(x.norm()))
					{
						cout << "M = " << endl << M << endl;
						cout << "f = " << endl << f << endl;
						assert(false);
					}*/
					for (HybridAlgebraicFace* ff : ce.RemovedFineFaces)
					{
						int localNumberInCE = ce.LocalRemovedFineFaceNumber(ff);
						for (int localCoarseFaceNumber = 0; localCoarseFaceNumber < ce.CoarseFaces.size(); ++localCoarseFaceNumber)
						{
							HybridFaceAggregate* cf = ce.CoarseFaces[localCoarseFaceNumber];
							chunk->Results.Coeffs.Add(ff->Number*fbs, cf->Number*fbs, x.block(localNumberInCE*fbs, localCoarseFaceNumber*fbs, fbs, fbs));
						}
					}
				});

			SparseMatrix InteriorThatReconstruct(Q_F->rows(), Q_F->cols());
			parallelLoop.Fill(InteriorThatReconstruct);


			NumberParallelLoop<CoeffsChunk> parallelLoop2(mesh.Faces.size());
			parallelLoop2.ReserveChunkCoeffsSize(fbs * 2 * fbs);
			parallelLoop2.Execute([this, fbs, &InteriorThatReconstruct, &Q_F, &mesh](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
				{
					const HybridAlgebraicFace* face = &mesh.Faces[faceNumber];
					if (face->IsRemovedOnCoarseMesh)
						chunk->Results.Coeffs.CopyRows(faceNumber*fbs, fbs, InteriorThatReconstruct);
					else
						chunk->Results.Coeffs.CopyRows(faceNumber*fbs, fbs, *Q_F);
				});
			P = new SparseMatrix(Q_F->rows(), Q_F->cols());
			parallelLoop2.Fill(*P);

		}
		else if (_multigridProlong == CAMGProlongation::HighOrder) // 8
		{
			SparseMatrix inv_A_T_Tc = Utils::InvertBlockDiagMatrix(*A_T_Tc, _cellBlockSize);
			SparseMatrix Theta = -inv_A_T_Tc * *A_T_Fc_tmp;   // Reconstruct
			SparseMatrix Pi = BuildHighOrderTraceOnRemovedFaces(mesh);
			
			SparseMatrix ReconstructAndTrace1 = Pi * Q_T * Theta;

			NumberParallelLoop<CoeffsChunk> parallelLoop(mesh.Faces.size());
			parallelLoop.ReserveChunkCoeffsSize(_faceBlockSize * 2 * _faceBlockSize);
			parallelLoop.Execute([this, &ReconstructAndTrace1, &Q_F, &mesh](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
				{
					const HybridAlgebraicFace* face = &mesh.Faces[faceNumber];
					if (face->IsRemovedOnCoarseMesh)
						chunk->Results.Coeffs.CopyRows(faceNumber*_faceBlockSize, _faceBlockSize, ReconstructAndTrace1);
					else
						chunk->Results.Coeffs.CopyRows(faceNumber*_faceBlockSize, _faceBlockSize, *Q_F);
				});
			P = new SparseMatrix(Q_F->rows(), Q_F->cols());
			parallelLoop.Fill(*P);
		}
		else
			Utils::FatalError("Unmanaged prolongation");

		//SparseMatrix* A_T_Fc = new SparseMatrix(Q_T.transpose() * A_T_F * *P); // Kills -prolong 1 or 2
		SparseMatrix* A_T_Fc = A_T_Fc_tmp;

		//SparseMatrix* A_F_Fc = new SparseMatrix(Q_F->transpose() * A_F_F * *Q_F);
		SparseMatrix* A_F_Fc = new SparseMatrix(P->transpose() * A_F_F.selfadjointView<Eigen::Lower>() * *P);


		//SparseMatrix* schurc = new SparseMatrix(ReduceSparsity(P->transpose() * schur * *P, mesh));
		SparseMatrix* schurc = new SparseMatrix(P->transpose() * schur.selfadjointView<Eigen::Lower>() * *P);

		return { A_T_Tc, A_T_Fc, A_F_Fc, P, Q_F, schurc, coarsestPossibleMeshReached };
	}

	void SetupDiscretizedOperator() override
	{
		/*SparseMatrix inv_A_T_T = Utils::InvertBlockDiagMatrix(*A_T_T, _cellBlockSize);
		SparseMatrix* schur = new SparseMatrix(*A_F_F - (A_T_F->transpose()) * (inv_A_T_T) * (*A_T_F));
		this->OperatorMatrix = schur;*/

		CondensedLevel* fine = dynamic_cast<CondensedLevel*>(this->FinerLevel);
		this->OperatorMatrix = new SparseMatrix(fine->Q_F.transpose() * *(fine->OperatorMatrix) * fine->Q_F);
	}

private:
	// Cell prolongation Q_T with only one 1 coefficient per row
	SparseMatrix BuildQ_T(const HybridAlgebraicMesh& mesh)
	{
		DenseMatrix Id = DenseMatrix::Identity(_cellBlockSize, _cellBlockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoopQ_T(mesh.Elements.size());
		parallelLoopQ_T.Execute([this, &mesh, &Id](BigNumber elemNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicElement& elem = mesh.Elements[elemNumber];
				chunk->Results.Coeffs.Add(elem.Number*_cellBlockSize, elem.CoarseElement->Number*_cellBlockSize, Id);
			});
		SparseMatrix Q_T = SparseMatrix(mesh.Elements.size()*_cellBlockSize, mesh.CoarseElements.size()*_cellBlockSize);
		parallelLoopQ_T.Fill(Q_T);
		return Q_T;
	}

	// Face prolongation FaceProlongation
	SparseMatrix BuildQ_F(const HybridAlgebraicMesh& mesh, const SparseMatrix& A_F_F)
	{
		bool enableAnisotropyManagement = false;
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);

		NumberParallelLoop<CoeffsChunk> parallelLoop1(mesh.Faces.size());
		parallelLoop1.Execute([this, &mesh, &Id, enableAnisotropyManagement, &A_F_F](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicFace* face = &mesh.Faces[faceNumber];
				if (face->IsRemovedOnCoarseMesh && (!Utils::ProgramArgs.Solver.MG.ManageAnisotropy || !face->CoarseFace))
				{
					// Take the average value of the coarse element faces
					HybridElementAggregate* elemAggreg = face->Elements[0]->CoarseElement;

					if (enableAnisotropyManagement)
					{
						map<HybridFaceAggregate*, double> couplings;
						double totalCouplings = 0;
						for (HybridFaceAggregate* coarseFace : elemAggreg->CoarseFaces)
						{
							double avgCouplingCoarseFace = 0;
							for (HybridAlgebraicFace* f : coarseFace->FineFaces)
							{
								DenseMatrix couplingBlock = A_F_F.block(face->Number*_faceBlockSize, f->Number*_faceBlockSize, _faceBlockSize, _faceBlockSize);
								double couplingFineFace = couplingBlock(0, 0);
								avgCouplingCoarseFace += couplingFineFace;
							}
							avgCouplingCoarseFace /= coarseFace->FineFaces.size();
							if (avgCouplingCoarseFace < 0)
							{
								couplings.insert({ coarseFace, avgCouplingCoarseFace });
								totalCouplings += avgCouplingCoarseFace;
							}
						}

						for (auto it = couplings.begin(); it != couplings.end(); it++)
						{
							HybridFaceAggregate* coarseFace = it->first;
							double coupling = it->second;
							chunk->Results.Coeffs.Add(face->Number, coarseFace->Number, -coupling / abs(totalCouplings) * Id);
						}
					}
					else
					{
						for (HybridFaceAggregate* coarseFace : elemAggreg->CoarseFaces)
							chunk->Results.Coeffs.Add(face->Number*_faceBlockSize, coarseFace->Number*_faceBlockSize, 1.0 / elemAggreg->CoarseFaces.size() * Id);
					}
				}
				else
					chunk->Results.Coeffs.Add(face->Number*_faceBlockSize, face->CoarseFace->Number*_faceBlockSize, Id);
			});


		SparseMatrix Q_F = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.CoarseFaces.size()*_faceBlockSize);
		parallelLoop1.Fill(Q_F);
		return Q_F;
	}

	SparseMatrix BuildQ_F_0Interior(const HybridAlgebraicMesh& mesh)
	{
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);

		NumberParallelLoop<CoeffsChunk> parallelLoop(mesh.Faces.size());
		parallelLoop.Execute([this, &mesh, &Id](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicFace* face = &mesh.Faces[faceNumber];
				if (!face->IsRemovedOnCoarseMesh)
					chunk->Results.Coeffs.Add(face->Number*_faceBlockSize, face->CoarseFace->Number*_faceBlockSize, Id);
			});


		SparseMatrix Q_F = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.CoarseFaces.size()*_faceBlockSize);
		parallelLoop.Fill(Q_F);
		return Q_F;
	}

	SparseMatrix BuildQ_F_AllAggregated(const AlgebraicMesh& skeleton)
	{
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoop(skeleton.Elements.size());
		parallelLoop.Execute([this, &skeleton, &Id](BigNumber elemNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const AlgebraicElement& elem = skeleton.Elements[elemNumber];
				chunk->Results.Coeffs.Add(elem.Number*_faceBlockSize, elem.CoarseElement->Number*_faceBlockSize, Id);
			});
		SparseMatrix Q_F = SparseMatrix(skeleton.Elements.size()*_faceBlockSize, skeleton.CoarseElements.size()*_faceBlockSize);
		parallelLoop.Fill(Q_F);
		return Q_F;
	}

	SparseMatrix BuildTrace(const HybridAlgebraicMesh& mesh)
	{
		// Pi: average on both sides of each face
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;

		NumberParallelLoop<CoeffsChunk> parallelLoopPi(mesh.Faces.size());
		parallelLoopPi.Execute([this, &mesh, &traceOfConstant](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicFace& face = mesh.Faces[faceNumber];
				if (face.IsRemovedOnCoarseMesh)
					chunk->Results.Coeffs.Add(faceNumber*_faceBlockSize, face.Elements[0]->Number*_cellBlockSize, traceOfConstant);
				else
				{
					assert(!face.Elements.empty());
					for (HybridAlgebraicElement* elem : face.Elements)
						chunk->Results.Coeffs.Add(faceNumber*_faceBlockSize, elem->Number*_cellBlockSize, 1.0 / face.Elements.size()*traceOfConstant);
				}
			});
		SparseMatrix Pi = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.Elements.size()*_cellBlockSize);
		parallelLoopPi.Fill(Pi);
		return Pi;
	}

	SparseMatrix BuildCoarseTraceOnFineRemovedFaces(const HybridAlgebraicMesh& mesh)
	{
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;

		NumberParallelLoop<CoeffsChunk> parallelLoopPi(mesh.Faces.size());
		parallelLoopPi.Execute([this, &mesh, &traceOfConstant](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicFace& face = mesh.Faces[faceNumber];
				if (face.IsRemovedOnCoarseMesh)
					chunk->Results.Coeffs.Add(faceNumber*_faceBlockSize, face.CoarseElements[0]->Number*_cellBlockSize, traceOfConstant);
			});
		SparseMatrix Pi = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.CoarseElements.size()*_cellBlockSize);
		parallelLoopPi.Fill(Pi);
		return Pi;
	}

	SparseMatrix BuildCoarseTraceOnFineFaces(const HybridAlgebraicMesh& mesh)
	{
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;

		NumberParallelLoop<CoeffsChunk> parallelLoopPi(mesh.Faces.size());
		parallelLoopPi.Execute([this, &mesh, &traceOfConstant](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicFace& face = mesh.Faces[faceNumber];
				for (HybridElementAggregate* ce : face.CoarseElements)
					chunk->Results.Coeffs.Add(faceNumber*_faceBlockSize, ce->Number*_cellBlockSize, (1.0 / face.CoarseElements.size())*traceOfConstant);
			});
		SparseMatrix Pi = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.CoarseElements.size()*_cellBlockSize);
		parallelLoopPi.Fill(Pi);
		return Pi;
	}

	SparseMatrix BuildHighOrderTraceOnRemovedFaces(const HybridAlgebraicMesh& mesh)
	{
		NumberParallelLoop<CoeffsChunk> parallelLoopPi(mesh.Faces.size());
		parallelLoopPi.Execute([this, &mesh](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicFace& face = mesh.Faces[faceNumber];
				if (face.IsRemovedOnCoarseMesh)
				{
					BigNumber elemNumber = face.Elements[0]->Number;
					DenseMatrix faceMass     = mesh.A_F_F->block(faceNumber * _faceBlockSize, faceNumber * _faceBlockSize, _faceBlockSize, _faceBlockSize);
					DenseMatrix cellFaceMass = mesh.A_T_F->block(elemNumber * _cellBlockSize, faceNumber * _faceBlockSize, _cellBlockSize, _faceBlockSize);
					DenseMatrix cellMass     = mesh.A_T_T->block(elemNumber * _cellBlockSize, elemNumber * _cellBlockSize, _cellBlockSize, _cellBlockSize);
					DenseMatrix trace = -faceMass.inverse() * cellFaceMass.transpose();
					//DenseMatrix trace = - cellFaceMass.transpose();
					chunk->Results.Coeffs.Add(faceNumber*_faceBlockSize, elemNumber*_cellBlockSize, trace);
				}
			});
		SparseMatrix Pi = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.Elements.size()*_cellBlockSize);
		parallelLoopPi.Fill(Pi);
		return Pi;
	}

	SparseMatrix ReduceSparsity(const SparseMatrix& A, const HybridAlgebraicMesh& mesh)
	{
		NumberParallelLoop<CoeffsChunk> parallelLoop(mesh.CoarseFaces.size());
		parallelLoop.Execute([this, &A, &mesh](BigNumber coarseFaceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridFaceAggregate& cf = mesh.CoarseFaces[coarseFaceNumber];

				for (int k = 0; k < _cellBlockSize; k++)
				{
					// RowMajor --> the following line iterates over the non-zeros of the elemNumber-th row.
					for (SparseMatrix::InnerIterator it(A, coarseFaceNumber*_faceBlockSize + k); it; ++it)
					{
						BigNumber coarseFaceNumber2 = it.col() / _faceBlockSize;
						const HybridFaceAggregate* cf2 = &mesh.CoarseFaces[coarseFaceNumber2];

						bool cf2IsInCf1Stencil = false;
						if (coarseFaceNumber == coarseFaceNumber2)
							cf2IsInCf1Stencil = true;
						else
						{
							for (HybridAlgebraicFace* ff : cf.FineFaces)
							{
								if (ff->IsRemovedOnCoarseMesh)
									continue;
								for (HybridElementAggregate* coarseElement : ff->CoarseElements)
								{
									if (find(coarseElement->CoarseFaces.begin(), coarseElement->CoarseFaces.end(), cf2) != coarseElement->CoarseFaces.end())
									{
										cf2IsInCf1Stencil = true;
										break;
									}
								}
								if (cf2IsInCf1Stencil)
									break;
							}
						}
						if (cf2IsInCf1Stencil)
							chunk->Results.Coeffs.Add(it.row(), it.col(), it.value());
					}
				}
			});

		SparseMatrix A2(A.rows(), A.cols());
		parallelLoop.Fill(A2);

		cout << "A.nonZeros()=" << A.nonZeros() << ", A2.nonZeros()=" << A2.nonZeros() << endl;

		return A2;
	}

public:
	void OnStartSetup() override
	{
		cout << "\t\tMesh                : " << this->A_T_T->rows() / _cellBlockSize << " elements, " << this->A_T_F->cols() / _faceBlockSize << " faces";
		if (!this->IsFinestLevel())
		{
			CondensedLevel* fine = dynamic_cast<CondensedLevel*>(this->FinerLevel);
			double nFine = fine->A_T_F->cols();
			double nCoarse = this->A_T_F->cols();
			cout << ", coarsening factor = " << (nFine/nCoarse);
		}
		cout << endl;
	}

	void SetupProlongation() override
	{}

	void SetupRestriction() override
	{
		double scalingFactor = 1.0;
		//scalingFactor = 1.0 / 4.0;
		R = (scalingFactor * P.transpose()).eval();
		//R = scalingFactor * P.transpose();
	}

	void OnEndSetup() override
	{
		if (this->CoarserLevel)
		{
			CondensedLevel* coarse = dynamic_cast<CondensedLevel*>(this->CoarserLevel);
			coarse->A_T_T = &A_T_Tc;
			coarse->A_T_F = &A_T_Fc;
			coarse->A_F_F = &A_F_Fc;
			coarse->inv_A_T_T = inv_A_T_Tc;
		}
	}

	~CondensedLevel()
	{
		if (!this->IsFinestLevel() && !this->UseGalerkinOperator)
			delete OperatorMatrix;
	}
};

class CondensedAMG : public Multigrid
{
private:
	CAMGFaceProlongation _faceProlong = CAMGFaceProlongation::FaceAggregates;
	CAMGProlongation _multigridProlong = CAMGProlongation::ReconstructionTrace1Step;
	int _cellBlockSize;
	int _faceBlockSize;
	double _strongCouplingThreshold;
public:

	CondensedAMG(int cellBlockSize, int faceBlockSize, double strongCouplingThreshold, CAMGFaceProlongation faceProlong, CAMGProlongation prolongation, int nLevels = 0)
		: Multigrid(nLevels)
	{
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
		this->_multigridProlong = prolongation;
		this->_faceProlong = faceProlong;
		this->BlockSizeForBlockSmoothers = faceBlockSize;
		this->UseGalerkinOperator = true;
		this->CoarseningStgy = CoarseningStrategy::DoublePairwiseAggregation;
		this->FaceCoarseningStgy = FaceCoarseningStrategy::InterfaceCollapsing;
		this->_fineLevel = new CondensedLevel(0, cellBlockSize, faceBlockSize, strongCouplingThreshold, faceProlong, prolongation);
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "CondensedAMG" << endl;

		os << "\t" << "Face prolongation       : ";
		if (_faceProlong == CAMGFaceProlongation::BoundaryAggregatesInteriorAverage)
			os << "aggregate interface faces, avg on interior faces ";
		else if (_faceProlong == CAMGFaceProlongation::BoundaryAggregatesInteriorZero)
			os << "aggregate interface faces, 0 on interior faces ";
		else if (_faceProlong == CAMGFaceProlongation::FaceAggregates)
			os << "aggregate all faces ";
		os << "[-face-prolong " << (unsigned)_faceProlong << "]" << endl;

		os << "\t" << "Prolongation            : ";
		if (_multigridProlong == CAMGProlongation::ReconstructionTrace1Step)
			os << "ReconstructionTrace1Step ";
		else if (_multigridProlong == CAMGProlongation::ReconstructionTrace2Steps)
			os << "ReconstructionTrace2Steps ";
		else if (_multigridProlong == CAMGProlongation::FaceProlongation)
			os << "FaceProlongation ";
		else if (_multigridProlong == CAMGProlongation::ReconstructTraceOrInject)
			os << "ReconstructTraceOrInject ";
		os << "[-prolong " << (unsigned)_multigridProlong << "]" << endl;
	}

	void EndSerialize(ostream& os) const override
	{
	}

	void Setup(const SparseMatrix& A) override
	{
		assert(false && "This Setup method cannot be used in this solver.");
	}

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F) override
	{
		CondensedLevel* fine = dynamic_cast<CondensedLevel*>(this->_fineLevel);
		fine->A_T_T = &A_T_T;
		fine->A_T_F = &A_T_F;
		fine->A_F_F = &A_F_F;
		SparseMatrix* inv_A_T_T = new SparseMatrix(Utils::InvertBlockDiagMatrix(A_T_T, _cellBlockSize));
		fine->inv_A_T_T = inv_A_T_T;
		Multigrid::Setup(A);
	}

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		if (initialGuessCode.compare("smooth") == 0)
			Utils::Warning("Smooth initial guess unmanaged in AMG.");
		return Multigrid::Solve(b, initialGuessCode);
	}
	
protected:
	Level* CreateCoarseLevel(Level* fineLevel) override
	{
		CondensedLevel* coarseLevel = new CondensedLevel(fineLevel->Number + 1, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold, _faceProlong, _multigridProlong);
		coarseLevel->OperatorMatrix = &dynamic_cast<CondensedLevel*>(fineLevel)->Ac;
		return coarseLevel;
	}
};