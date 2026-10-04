#pragma once
#include <mutex>
#include <unsupported/Eigen/SparseExtra>
#include "HybridAlgebraicMesh.h"
#include "../AggregAMG/AlgebraicMesh.h"
#include "../Level.h"
#include "../../../Utils/SparseMatrixOps.h"
using namespace std;

class UncondensedLevel : public Level
{
private:
	UAMGFaceProlongation _faceProlong = UAMGFaceProlongation::FaceAggregates;
	UAMGProlongation _coarseningProlong = UAMGProlongation::FaceProlongation;
	UAMGProlongation _multigridProlong = UAMGProlongation::ReconstructSmoothedTraceOrInject;
	int _degree;
	int _cellBlockSize;
	int _faceBlockSize;
	double _strongCouplingThreshold;
public:
	const SparseMatrix* A_T_T;
	const SparseMatrix* A_T_F;
	const SparseMatrix* A_F_F;
private:
	//SparseMatrix Q_T;
	SparseMatrix Q_F;

	SparseMatrix A_T_Tc;
	SparseMatrix A_T_Fc;
	SparseMatrix A_F_Fc;
public:
	SparseMatrix Ac;

	// Set by UncondensedAMG, to skip the computations that the chosen options don't use:
	// - the coarse blocks A_F_F, computed at each coarsening step;
	// - the face prolongation Q_F chained over the coarsening steps.
	// If not computed, the coarse level's A_F_F is null.
	bool ComputeCoarseA_F_F = true;
	bool ComputeQ_F = true;

public:
	UncondensedLevel(int number, int degree, int cellBlockSize, int faceBlockSize, double strongCouplingThreshold, UAMGFaceProlongation faceProlong, UAMGProlongation coarseningProlong, UAMGProlongation mgProlong)
		: Level(number)
	{
		//if (cellBlockSize > 1 || faceBlockSize > 1)
			//Utils::Warning("This multigrid is efficient if cellBlockSize = faceBlockSize = 1. It may converge badly.");

		this->_degree = degree;
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
		this->_faceProlong = faceProlong;
		this->_coarseningProlong = coarseningProlong;
		this->_multigridProlong = mgProlong;
	}

	BigNumber NUnknowns() override
	{
		return A_T_F->cols();
	}

	int PolynomialDegree() override
	{
		return _degree;
	}

	int BlockSizeForBlockSmoothers() override
	{
		return _faceBlockSize;
	}

	int CellBlockSize()
	{
		return _cellBlockSize;
	}

	int FaceBlockSize()
	{
		return _faceBlockSize;
	}

	void CoarsenMesh(H_CoarsStgy coarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, FaceCollapsing bdryFaceCollapsing, double requestedCoarseningFactor, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) override
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

		SparseMatrix *auxP, *auxQ_F, *auxSchur;
		HybridAlgebraicMesh *mesh, *coarseMesh;
		const SparseMatrix* schur = this->OperatorMatrix;
		double actualCoarseningFactor = 0;
		double nCoarsenings = 0;

		HybridAlgebraicMesh initialFineMesh(A_T_T, A_T_F, A_F_F, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold);
		mesh = &initialFineMesh;

		while (!CoarseningCriteriaReached(coarseningStgy, requestedCoarseningFactor, nCoarsenings, actualCoarseningFactor))
		{
			// Coarsening
			std::tie(coarseMesh, auxP, auxQ_F, auxSchur, coarsestPossibleMeshReached) = Coarsen(*mesh, *schur, coarseningStgy, faceCoarseningStgy);
			if (coarsestPossibleMeshReached)
				return;

			// Global prolongation operator
			if (nCoarsenings == 0)
			{
				if (_multigridProlong == UAMGProlongation::ChainedCoarseningProlongations)
					this->P = *auxP;
				if (ComputeQ_F)
					this->Q_F = *auxQ_F;
			}
			else
			{
				if (_multigridProlong == UAMGProlongation::ChainedCoarseningProlongations)
					this->P = SparseMatrixOps::Multiply(this->P, *auxP);
				if (ComputeQ_F)
					this->Q_F = SparseMatrixOps::Multiply(this->Q_F, *auxQ_F);
			}
			delete auxP;
			if (auxQ_F != auxP)
				delete auxQ_F;

			// Global coarsening factor
			double nFine = initialFineMesh.A_T_F->cols();
			double nCoarse = coarseMesh->A_T_F->cols();
			actualCoarseningFactor = nFine / nCoarse;
			cout << "\t\tCoarsening factor = " << actualCoarseningFactor << endl;
			if (actualCoarseningFactor >= 10)
				Utils::Warning("The coarsening seems a little too strong...");


			// Prepare next coarsening
			if (nCoarsenings > 0)
			{
				if (_multigridProlong != UAMGProlongation::ChainedCoarseningProlongations)
				{
					// Update global coarsening
					UpdateGlobalCoarsening(initialFineMesh, *mesh);
				}

				// Print aggregates for debugging purposes
				/*UpdateGlobalCoarsening(initialFineMesh, *mesh);
				for (HybridElementAggregate& ce : initialFineMesh.CoarseElements)
				{
					cout << "Aggregate " << ce.Number << ": ";
					for (HybridAlgebraicElement* fe : ce.FineElements)
						cout << fe->Number << " ";
					cout << endl;
				}*/

				// Memory release
				delete mesh->A_T_T;
				delete mesh->A_T_F;
				if (mesh->A_F_F) delete mesh->A_F_F;
				delete schur;
				delete mesh;
			}

			mesh = coarseMesh;
			schur = auxSchur;

			nCoarsenings++;
		}

		// Ending...
		if (_multigridProlong == UAMGProlongation::ChainedCoarseningProlongations)
		{
			this->A_T_Tc = std::move(*coarseMesh->A_T_T);
			this->A_T_Fc = std::move(*coarseMesh->A_T_F);
			if (coarseMesh->A_F_F)
				this->A_F_Fc = std::move(*coarseMesh->A_F_F);
			this->Ac     = std::move(*schur);
		}
		else
		{
			// Multigrid prolongation
			SparseMatrix Q_T = BuildQ_T(initialFineMesh);
			coarseMesh->BuildElementFaces(); // the coarse cells and faces are used by the reconstruction (Theta)
			SparseMatrix* P = BuildProlongation(_multigridProlong, initialFineMesh, *this->OperatorMatrix, *coarseMesh, Q_T, &Q_F);
			this->P = std::move(*P);
			delete P;

			this->A_T_Tc = std::move(*coarseMesh->A_T_T);
			SparseMatrix Pt = SparseMatrixOps::Transpose(this->P);
			this->A_T_Fc = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(SparseMatrixOps::Transpose(Q_T), *this->A_T_F), this->P);
			if (ComputeCoarseA_F_F)
				this->A_F_Fc = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, *this->A_F_F), this->P);
			this->Ac     = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, *this->OperatorMatrix), this->P);
		}
		delete coarseMesh;
		delete schur;
	}

	void UpdateGlobalCoarsening(HybridAlgebraicMesh& initialFineMesh, HybridAlgebraicMesh& currentMesh)
	{
		#pragma omp parallel for
		for (BigNumber aggregNumber = 0; aggregNumber < currentMesh.CoarseElements.size(); ++aggregNumber)
		{
			HybridElementAggregate& aggreg = currentMesh.CoarseElements[aggregNumber];
			RemoveInitialFineFaces(initialFineMesh, aggreg);
			UpdateInitialFineElements(initialFineMesh, aggreg);
		}

		#pragma omp parallel for
		for (BigNumber aggregNumber = 0; aggregNumber < currentMesh.CoarseFaces.size(); ++aggregNumber)
		{
			HybridFaceAggregate& aggreg = currentMesh.CoarseFaces[aggregNumber];
			UpdateRemainingInitialFineFaces(initialFineMesh, aggreg);
		}

		initialFineMesh.CoarseElements = std::move(currentMesh.CoarseElements);
		initialFineMesh.CoarseFaces = std::move(currentMesh.CoarseFaces);
	}

	void UpdateInitialFineElements(HybridAlgebraicMesh& initialFineMesh, HybridElementAggregate& aggreg)
	{
		vector<HybridAlgebraicElement*> fineElements;
		for (HybridAlgebraicElement* ce : aggreg.FineElements)
		{
			HybridElementAggregate& finerAggreg = initialFineMesh.CoarseElements[ce->Number];
			for (HybridAlgebraicElement* fe : finerAggreg.FineElements)
			{
				fineElements.push_back(fe);
				fe->CoarseElement = &aggreg;
			}
		}
		aggreg.FineElements = fineElements;
	}

	void UpdateRemainingInitialFineFaces(HybridAlgebraicMesh& initialFineMesh, HybridFaceAggregate& aggreg)
	{
		vector<HybridAlgebraicFace*> fineFaces;
		for (HybridAlgebraicFace* cf : aggreg.FineFaces)
		{
			HybridFaceAggregate& finerAggreg = initialFineMesh.CoarseFaces[cf->Number];
			for (HybridAlgebraicFace* ff : finerAggreg.FineFaces)
			{
				assert(!ff->IsRemovedOnCoarseMesh);
				fineFaces.push_back(ff);
				ff->CoarseFace = &aggreg;
				ff->CoarseElements = {};
			}
		}
		aggreg.FineFaces = fineFaces;
	}

	void RemoveInitialFineFaces(HybridAlgebraicMesh& initialFineMesh, HybridElementAggregate& aggreg)
	{
		vector<HybridAlgebraicFace*> removedFineFaces;
		for (HybridAlgebraicFace* cf : aggreg.RemovedFineFaces)
		{
			HybridFaceAggregate& finerAggreg = initialFineMesh.CoarseFaces[cf->Number];
			for (HybridAlgebraicFace* ff : finerAggreg.FineFaces)
			{
				assert(!ff->IsRemovedOnCoarseMesh);
				ff->IsRemovedOnCoarseMesh = true;
				ff->CoarseFace = nullptr;
				ff->CoarseElements = { &aggreg };
				removedFineFaces.push_back(ff);
			}
		}
		for (HybridAlgebraicElement* ce : aggreg.FineElements)
		{
			HybridElementAggregate& finerAggreg = initialFineMesh.CoarseElements[ce->Number];
			for (HybridAlgebraicFace* ff : finerAggreg.RemovedFineFaces)
			{
				assert(ff->IsRemovedOnCoarseMesh);
				assert(!ff->CoarseFace);
				ff->CoarseElements = { &aggreg };
				removedFineFaces.push_back(ff);
			}
		}
		aggreg.RemovedFineFaces = removedFineFaces;
	}


	bool CoarseningCriteriaReached(H_CoarsStgy coarseningStgy, double requestedCoarseningRatio, int nCoarseningsPerformed, double coarseningRatio)
	{
		if (coarseningStgy == H_CoarsStgy::DoublePairwiseAggregation)
			return nCoarseningsPerformed == 2;
		if (coarseningStgy == H_CoarsStgy::MultiplePairwiseAggregation)
			return coarseningRatio >= requestedCoarseningRatio;
		if (coarseningStgy == H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours)
			return nCoarseningsPerformed == 1; // only 1 pass of coarsening
		if (coarseningStgy == H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours)
			return coarseningRatio >= requestedCoarseningRatio;
		return nCoarseningsPerformed == 1;
	}




	// Returns <coarseMesh, P, Q_F, schurc, coarsestPossibleMeshReached>
	tuple<HybridAlgebraicMesh*, SparseMatrix*, SparseMatrix*, SparseMatrix*, bool> Coarsen(HybridAlgebraicMesh& mesh, const SparseMatrix& schur, H_CoarsStgy elemCoarseningStgy, FaceCoarseningStrategy faceCoarseningStgy)
	{
		//ExportMatrix(A_T_T, "A_T_T", 0);
		//ExportMatrix(A_T_F, "A_T_F", 0);
		//ExportMatrix(A_F_F, "A_F_F", 0);

		bool coarsestPossibleMeshReached = false;

		bool onlyFacesUsed = this->_faceProlong == UAMGFaceProlongation::FaceAggregates && _multigridProlong == UAMGProlongation::FaceProlongation;

		if (!onlyFacesUsed)
		{
			mesh.Build();
			mesh.Coarsen(elemCoarseningStgy, faceCoarseningStgy, coarsestPossibleMeshReached);
			if (coarsestPossibleMeshReached)
				return { nullptr, nullptr, nullptr, nullptr, coarsestPossibleMeshReached };
		}

		// Cell-prolongation operator with only one 1 coefficient per row
		SparseMatrix Q_T = BuildQ_T(mesh);

		// Face-prolongation operator
		SparseMatrix* Q_F;
		if (this->_faceProlong == UAMGFaceProlongation::BoundaryAggregatesInteriorAverage)
		{
			// Face-prolongation operator with only one 1 coefficient per row for kept or aggregated faces, average for removed faces
			Q_F = new SparseMatrix(BuildQ_F(mesh));
		}
		else if (this->_faceProlong == UAMGFaceProlongation::BoundaryAggregatesInteriorZero)
		{
			Q_F = new SparseMatrix(BuildQ_F_0Interior(mesh));
		}
		else if (this->_faceProlong == UAMGFaceProlongation::FaceAggregates)
		{
			AlgebraicMesh skeleton(_faceBlockSize, 0);
			//skeleton.Build(*A_F_F);
			skeleton.Build(schur);
			skeleton.PairWiseAggregate(coarsestPossibleMeshReached);
			if (coarsestPossibleMeshReached)
				return { nullptr, nullptr, nullptr, nullptr, coarsestPossibleMeshReached };
			Q_F = new SparseMatrix(BuildQ_F_AllAggregated(skeleton));
		}
		else
			Utils::FatalError("Unmanaged -face-prolong");

		// Intermediate coarse operators
		// The products are evaluated from left to right, and the symmetric matrices are given by their lower triangular part
		SparseMatrix Q_Tt = SparseMatrixOps::Transpose(Q_T);
		SparseMatrix Q_Tt_A_T_F = SparseMatrixOps::Multiply(Q_Tt, *mesh.A_T_F);
		SparseMatrix* A_T_Tc = nullptr;
		SparseMatrix* A_T_Fc_tmp = nullptr;
		if (!onlyFacesUsed)
		{
			A_T_Tc     = new SparseMatrix(SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Q_Tt, SparseMatrixOps::FullFromLower(*mesh.A_T_T)), Q_T));
			A_T_Fc_tmp = new SparseMatrix(SparseMatrixOps::Multiply(Q_Tt_A_T_F, *Q_F));
		}

		// The prolongation only uses its matrices and the faces of its elements (reconstruction, Theta())
		HybridAlgebraicMesh auxCoarseMesh(A_T_Tc, A_T_Fc_tmp, nullptr, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold);
		auxCoarseMesh.BuildElementFaces();

		/*ExportMatrix(Q_T, "Q_T", 0);
		ExportMatrix(*Q_F, "Q_F", 0);
		ExportMatrix(*A_T_Tc, "A_T_Tc", 0);
		ExportMatrix(*A_T_Fc_tmp, "A_T_Fc", 0);*/

		// Multigrid prolongation
		SparseMatrix* P = BuildProlongation(this->_coarseningProlong, mesh, schur, auxCoarseMesh, Q_T, Q_F);

		SparseMatrix* A_T_Fc = new SparseMatrix(SparseMatrixOps::Multiply(Q_Tt_A_T_F, *P)); // Kills -prolong 1 or 2 because P is then very dense
		//SparseMatrix* A_T_Fc = A_T_Fc_tmp;

		SparseMatrix Pt = SparseMatrixOps::Transpose(*P);
		SparseMatrix* A_F_Fc = nullptr;
		if (ComputeCoarseA_F_F)
			A_F_Fc = new SparseMatrix(SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, SparseMatrixOps::FullFromLower(*mesh.A_F_F)), *P));
		SparseMatrix* schurc = new SparseMatrix(SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, SparseMatrixOps::FullFromLower(schur)), *P));

		HybridAlgebraicMesh* coarseMesh = new HybridAlgebraicMesh(A_T_Tc, A_T_Fc, A_F_Fc, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold);
		
		return { coarseMesh, P, Q_F, schurc, coarsestPossibleMeshReached };
	}



	SparseMatrix* BuildProlongation(UAMGProlongation prolong, HybridAlgebraicMesh& mesh, const SparseMatrix& schur, HybridAlgebraicMesh& coarseMesh,
									const SparseMatrix& Q_T, SparseMatrix* Q_F)
	{
		SparseMatrix* P;
		if (prolong == UAMGProlongation::ReconstructionTrace) // 1
		{
			// Theta: reconstruction from the coarse faces to the coarse cells
			SparseMatrix Theta = coarseMesh.Theta();
			// Pi: average on both sides of each face
			SparseMatrix Pi = BuildTrace(mesh);

			P = new SparseMatrix(SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pi, Q_T), Theta));
		}
		else if (prolong == UAMGProlongation::FaceProlongation) // 3
		{
			// -g 1 -prolong 3 -face-prolong 3 -cs z
			P = Q_F;
		}
		else if (prolong == UAMGProlongation::FaceProlongationAndInteriorSmoothing) // 4
		{
			vector<bool> isRemoved = RemovedFaces(mesh);

			// Smoothing (only the rows of the removed faces are kept)
			BlockJacobi blockJacobi(_faceBlockSize, 2.0 / 3.0);
			blockJacobi.Setup(schur);
			SparseMatrix J = blockJacobi.IterationMatrix(isRemoved);

			SparseMatrix smoothedQ_F = SparseMatrixOps::Multiply(J, *Q_F);

			P = new SparseMatrix(RemovedOrKeptFaceRows(isRemoved, smoothedQ_F, *Q_F));
		}
		else if (prolong == UAMGProlongation::ReconstructTraceOrInject) // 5
		{
			SparseMatrix Theta = coarseMesh.Theta();   // Reconstruct
			SparseMatrix Pi = BuildCoarseTraceOnFineRemovedFaces(mesh); // Trace

			SparseMatrix ReconstructAndTrace = SparseMatrixOps::Multiply(Pi, Theta);

			P = new SparseMatrix(RemovedOrKeptFaceRows(RemovedFaces(mesh), ReconstructAndTrace, *Q_F));
		}
		else if (prolong == UAMGProlongation::ReconstructSmoothedTraceOrInject) // 6
		{
			vector<bool> isRemoved = RemovedFaces(mesh);

			SparseMatrix Theta = coarseMesh.Theta();   // Reconstruct
			SparseMatrix Pi = BuildCoarseTraceOnFineRemovedFaces(mesh); // Trace

			SparseMatrix ReconstructAndTrace = SparseMatrixOps::Multiply(Pi, Theta);

			SparseMatrix ReconstructTraceOrInject = RemovedOrKeptFaceRows(isRemoved, ReconstructAndTrace, *Q_F);

			// Smoothing (only the rows of the removed faces are kept: the matrix J is not assembled for the others)
			BlockJacobi blockJacobi(_faceBlockSize, 2.0/3.0);
			blockJacobi.Setup(schur);
			SparseMatrix J = blockJacobi.IterationMatrix(isRemoved);

			SparseMatrix ReconstructAndSmoothedTrace = SparseMatrixOps::Multiply(J, ReconstructTraceOrInject);

			P = new SparseMatrix(RemovedOrKeptFaceRows(isRemoved, ReconstructAndSmoothedTrace, *Q_F));
		}
		else if (prolong == UAMGProlongation::FindInteriorThatReconstructs) // 7
		{
			int cbs = _cellBlockSize;
			int fbs = _faceBlockSize;

			ThreadLocalCoeffs coeffs;
			#pragma omp parallel for
			for (BigNumber ceNumber = 0; ceNumber < mesh.CoarseElements.size(); ++ceNumber)
			{
				const HybridElementAggregate& ce = mesh.CoarseElements[ceNumber];

				if (ce.RemovedFineFaces.empty())
					continue;

				// Construction of Theta_Tc
				DenseMatrix A_Tc_Tc = coarseMesh.A_T_T->block(ce.Number, ce.Number, cbs, cbs);

				DenseMatrix A_Tc_F(cbs, ce.CoarseFaces.size()*fbs);
				for (int cfLocalNumber = 0; cfLocalNumber < ce.CoarseFaces.size(); cfLocalNumber++)
				{
					HybridFaceAggregate* cf = ce.CoarseFaces[cfLocalNumber];
					A_Tc_F.block(0, cfLocalNumber*fbs, cbs, fbs) = coarseMesh.A_T_F->block(ce.Number*cbs, cf->Number*fbs, cbs, fbs);
				}
				DenseMatrix Theta_Tc = -A_Tc_Tc.llt().solve(A_Tc_F);

				// Matrix for minimization problem
				DenseMatrix M(ce.FineElements.size()*cbs, ce.RemovedFineFaces.size()*fbs);

				// RHS for minimization problem
				DenseMatrix f = DenseMatrix::Zero(ce.FineElements.size()*cbs, ce.CoarseFaces.size()*fbs);

				for (int feLocalNumber = 0; feLocalNumber < ce.FineElements.size(); ++feLocalNumber)
				{
					HybridAlgebraicElement* fe = ce.FineElements[feLocalNumber];

					// Construction of Theta_Tf
					DenseMatrix A_Tf_Tf = mesh.A_T_T->block(fe->Number, fe->Number, cbs, cbs);
					DenseMatrix A_Tf_F(cbs, fe->Faces.size()*fbs);
					for (int ffLocalNumber = 0; ffLocalNumber < fe->Faces.size(); ffLocalNumber++)
					{
						HybridAlgebraicFace* ff = fe->Faces[ffLocalNumber];
						A_Tf_F.block(0, ffLocalNumber*fbs, cbs, fbs) = mesh.A_T_F->block(fe->Number*cbs, ff->Number*fbs, cbs, fbs);
					}
					DenseMatrix Theta_Tf = -A_Tf_Tf.llt().solve(A_Tf_F);
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
						coeffs.Local().Add(ff->Number*fbs, cf->Number*fbs, x.block(localNumberInCE*fbs, localCoarseFaceNumber*fbs, fbs, fbs));
					}
				}
			}

			SparseMatrix InteriorThatReconstruct(Q_F->rows(), Q_F->cols());
			coeffs.Fill(InteriorThatReconstruct);

			P = new SparseMatrix(RemovedOrKeptFaceRows(RemovedFaces(mesh), InteriorThatReconstruct, *Q_F));
		}
		else if (prolong == UAMGProlongation::HighOrder) // 8
		{
			SparseMatrix Theta = coarseMesh.Theta();
			SparseMatrix Pi = BuildHighOrderTraceOnRemovedFaces(mesh);

			SparseMatrix ReconstructAndTrace1 = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pi, Q_T), Theta);

			P = new SparseMatrix(RemovedOrKeptFaceRows(RemovedFaces(mesh), ReconstructAndTrace1, *Q_F));
		}
		else if (prolong == UAMGProlongation::ReconstructionTranspose2Steps) // 9
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
		else
			Utils::FatalError("Unmanaged prolongation");

		return P;
	}

	void SetupDiscretizedOperator() override
	{
		/*SparseMatrix inv_A_T_T = Utils::InvertBlockDiagMatrix(*A_T_T, _cellBlockSize);
		SparseMatrix* schur = new SparseMatrix(*A_F_F - (A_T_F->transpose()) * (inv_A_T_T) * (*A_T_F));
		this->OperatorMatrix = schur;*/

		UncondensedLevel* fine = dynamic_cast<UncondensedLevel*>(this->FinerLevel);
		this->OperatorMatrix = new SparseMatrix(SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(SparseMatrixOps::Transpose(fine->Q_F), *(fine->OperatorMatrix)), fine->Q_F));
	}

public:
	void SetupOperatorByBlockExtraction()
	{
		UncondensedLevel* fine = dynamic_cast<UncondensedLevel*>(FinerLevel);
		if (!fine->A_F_F)
			Utils::FatalError("UncondensedAMG: the p-coarsening requires the block A_F_F, which has not been computed on the finer level (see UncondensedAMG::CoarseA_F_FNeeded()).");

		auto nElems = fine->A_T_F->rows() / fine->_cellBlockSize;
		auto nFaces = fine->A_T_F->cols() / fine->_faceBlockSize;

		this->OperatorMatrix = ExtractCoarseMatrix(*this->FinerLevel->OperatorMatrix, nFaces, nFaces, fine->_faceBlockSize, fine->_faceBlockSize, this->_faceBlockSize, this->_faceBlockSize);
		this->A_T_T = ExtractCoarseMatrix(*fine->A_T_T, nElems, nElems, fine->_cellBlockSize, fine->_cellBlockSize, this->_cellBlockSize, this->_cellBlockSize);
		this->A_T_F = ExtractCoarseMatrix(*fine->A_T_F, nElems, nFaces, fine->_cellBlockSize, fine->_faceBlockSize, this->_cellBlockSize, this->_faceBlockSize);
		this->A_F_F = ExtractCoarseMatrix(*fine->A_F_F, nFaces, nFaces, fine->_faceBlockSize, fine->_faceBlockSize, this->_faceBlockSize, this->_faceBlockSize);
	}

private:
	static SparseMatrix* ExtractCoarseMatrix(const SparseMatrix& fineMatrix, BigNumber nBlockRows, BigNumber nBlockCols, int fineBlockRows, int fineBlockCols, int coarseBlockRows, int coarseBlockCols)
	{
		ThreadLocalCoeffs coeffs;
		#pragma omp parallel for
		for (BigNumber i = 0; i < nBlockRows; ++i)
		{
			for (int k = 0; k < coarseBlockRows; k++)
			{
				for (RowMajorSparseMatrix::InnerIterator it(fineMatrix, i*fineBlockRows + k); it; ++it)
				{
					auto j = it.col() / fineBlockCols;
					int l = it.col() - j * fineBlockCols;
					if (l < coarseBlockCols)
						coeffs.Local().Add(i*coarseBlockRows + k, j*coarseBlockCols + l, it.value());
				}
			}
		}
		SparseMatrix* coarseMatrix = new SparseMatrix(nBlockRows * coarseBlockRows, nBlockCols * coarseBlockCols);
		coeffs.Fill(*coarseMatrix);
		return coarseMatrix;
	}


	// For each face of the mesh, whether it is removed on the coarse mesh (i.e. interior to an aggregate)
	static vector<bool> RemovedFaces(const HybridAlgebraicMesh& mesh)
	{
		vector<bool> isRemoved(mesh.Faces.size());
		for (BigNumber faceNumber = 0; faceNumber < mesh.Faces.size(); ++faceNumber)
			isRemoved[faceNumber] = mesh.Faces[faceNumber].IsRemovedOnCoarseMesh;
		return isRemoved;
	}

	// Prolongation made of the rows of removedFaceRows for the removed faces, and of keptFaceRows for the others.
	// The coefficients of absolute value <= NonZeroCoefficients::ZeroThreshold are dropped.
	SparseMatrix RemovedOrKeptFaceRows(const vector<bool>& isRemoved, const SparseMatrix& removedFaceRows, const SparseMatrix& keptFaceRows)
	{
		return SparseMatrixOps::SelectRows(isRemoved, removedFaceRows, keptFaceRows, _faceBlockSize, NonZeroCoefficients::ZeroThreshold);
	}

	// Cell prolongation Q_T with only one 1 coefficient per row
	SparseMatrix BuildQ_T(const HybridAlgebraicMesh& mesh)
	{
		DenseMatrix Id = DenseMatrix::Identity(_cellBlockSize, _cellBlockSize);
		ThreadLocalCoeffs coeffsQ_T;
		#pragma omp parallel for
		for (BigNumber elemNumber = 0; elemNumber < mesh.Elements.size(); ++elemNumber)
		{
			const HybridAlgebraicElement& elem = mesh.Elements[elemNumber];
			coeffsQ_T.Local().Add(elem.Number*_cellBlockSize, elem.CoarseElement->Number*_cellBlockSize, Id);
		}
		SparseMatrix Q_T = SparseMatrix(mesh.Elements.size()*_cellBlockSize, mesh.CoarseElements.size()*_cellBlockSize);
		coeffsQ_T.Fill(Q_T);
		return Q_T;
	}

	// Face prolongation FaceProlongation
	SparseMatrix BuildQ_F(const HybridAlgebraicMesh& mesh)
	{
		bool enableAnisotropyManagement = false;
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);

		ThreadLocalCoeffs coeffs1;
		#pragma omp parallel for
		for (BigNumber faceNumber = 0; faceNumber < mesh.Faces.size(); ++faceNumber)
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
							DenseMatrix couplingBlock = mesh.A_F_F->block(face->Number*_faceBlockSize, f->Number*_faceBlockSize, _faceBlockSize, _faceBlockSize);
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
						coeffs1.Local().Add(face->Number, coarseFace->Number, -coupling / abs(totalCouplings) * Id);
					}
				}
				else
				{
					for (HybridFaceAggregate* coarseFace : elemAggreg->CoarseFaces)
						coeffs1.Local().Add(face->Number*_faceBlockSize, coarseFace->Number*_faceBlockSize, 1.0 / elemAggreg->CoarseFaces.size() * Id);
				}
			}
			else
				coeffs1.Local().Add(face->Number*_faceBlockSize, face->CoarseFace->Number*_faceBlockSize, Id);
		}


		SparseMatrix Q_F = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.CoarseFaces.size()*_faceBlockSize);
		coeffs1.Fill(Q_F);
		return Q_F;
	}

	SparseMatrix BuildQ_F_0Interior(const HybridAlgebraicMesh& mesh)
	{
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);

		ThreadLocalCoeffs coeffs;
		#pragma omp parallel for
		for (BigNumber faceNumber = 0; faceNumber < mesh.Faces.size(); ++faceNumber)
		{
			const HybridAlgebraicFace* face = &mesh.Faces[faceNumber];
			if (!face->IsRemovedOnCoarseMesh)
				coeffs.Local().Add(face->Number*_faceBlockSize, face->CoarseFace->Number*_faceBlockSize, Id);
		}


		SparseMatrix Q_F = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.CoarseFaces.size()*_faceBlockSize);
		coeffs.Fill(Q_F);
		return Q_F;
	}

	SparseMatrix BuildQ_F_AllAggregated(const AlgebraicMesh& skeleton)
	{
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);
		ThreadLocalCoeffs coeffs;
		#pragma omp parallel for
		for (BigNumber elemNumber = 0; elemNumber < skeleton.Elements.size(); ++elemNumber)
		{
			const AlgebraicElement& elem = skeleton.Elements[elemNumber];
			coeffs.Local().Add(elem.Number*_faceBlockSize, elem.CoarseElement->Number*_faceBlockSize, Id);
		}
		SparseMatrix Q_F = SparseMatrix(skeleton.Elements.size()*_faceBlockSize, skeleton.CoarseElements.size()*_faceBlockSize);
		coeffs.Fill(Q_F);
		return Q_F;
	}

	SparseMatrix BuildTrace(const HybridAlgebraicMesh& mesh)
	{
		// Pi: average on both sides of each face
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;

		ThreadLocalCoeffs coeffsPi;
		#pragma omp parallel for
		for (BigNumber faceNumber = 0; faceNumber < mesh.Faces.size(); ++faceNumber)
		{
			const HybridAlgebraicFace& face = mesh.Faces[faceNumber];
			if (face.IsRemovedOnCoarseMesh)
				coeffsPi.Local().Add(faceNumber*_faceBlockSize, face.Elements[0]->Number*_cellBlockSize, traceOfConstant);
			else
			{
				assert(!face.Elements.empty());
				for (HybridAlgebraicElement* elem : face.Elements)
					coeffsPi.Local().Add(faceNumber*_faceBlockSize, elem->Number*_cellBlockSize, 1.0 / face.Elements.size()*traceOfConstant);
			}
		}
		SparseMatrix Pi = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.Elements.size()*_cellBlockSize);
		coeffsPi.Fill(Pi);
		return Pi;
	}

	SparseMatrix BuildCoarseTraceOnFineRemovedFaces(const HybridAlgebraicMesh& mesh)
	{
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;

		ThreadLocalCoeffs coeffsPi;
		#pragma omp parallel for
		for (BigNumber faceNumber = 0; faceNumber < mesh.Faces.size(); ++faceNumber)
		{
			const HybridAlgebraicFace& face = mesh.Faces[faceNumber];
			if (face.IsRemovedOnCoarseMesh)
				coeffsPi.Local().Add(faceNumber*_faceBlockSize, face.CoarseElements[0]->Number*_cellBlockSize, traceOfConstant);
		}
		SparseMatrix Pi = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.CoarseElements.size()*_cellBlockSize);
		coeffsPi.Fill(Pi);
		return Pi;
	}

	SparseMatrix BuildCoarseTraceOnFineFaces(const HybridAlgebraicMesh& mesh)
	{
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;

		ThreadLocalCoeffs coeffsPi;
		#pragma omp parallel for
		for (BigNumber faceNumber = 0; faceNumber < mesh.Faces.size(); ++faceNumber)
		{
			const HybridAlgebraicFace& face = mesh.Faces[faceNumber];
			for (HybridElementAggregate* ce : face.CoarseElements)
				coeffsPi.Local().Add(faceNumber*_faceBlockSize, ce->Number*_cellBlockSize, (1.0 / face.CoarseElements.size())*traceOfConstant);
		}
		SparseMatrix Pi = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.CoarseElements.size()*_cellBlockSize);
		coeffsPi.Fill(Pi);
		return Pi;
	}

	SparseMatrix BuildHighOrderTraceOnRemovedFaces(const HybridAlgebraicMesh& mesh)
	{
		ThreadLocalCoeffs coeffsPi;
		#pragma omp parallel for
		for (BigNumber faceNumber = 0; faceNumber < mesh.Faces.size(); ++faceNumber)
		{
			const HybridAlgebraicFace& face = mesh.Faces[faceNumber];
			if (face.IsRemovedOnCoarseMesh)
			{
				BigNumber elemNumber = face.Elements[0]->Number;
				DenseMatrix faceMass     = mesh.A_F_F->block(faceNumber * _faceBlockSize, faceNumber * _faceBlockSize, _faceBlockSize, _faceBlockSize);
				DenseMatrix cellFaceMass = mesh.A_T_F->block(elemNumber * _cellBlockSize, faceNumber * _faceBlockSize, _cellBlockSize, _faceBlockSize);
				DenseMatrix cellMass     = mesh.A_T_T->block(elemNumber * _cellBlockSize, elemNumber * _cellBlockSize, _cellBlockSize, _cellBlockSize);
				DenseMatrix trace = -faceMass.llt().solve(cellFaceMass.transpose());
				//DenseMatrix trace = - cellFaceMass.transpose();
				coeffsPi.Local().Add(faceNumber*_faceBlockSize, elemNumber*_cellBlockSize, trace);
			}
		}
		SparseMatrix Pi = SparseMatrix(mesh.Faces.size()*_faceBlockSize, mesh.Elements.size()*_cellBlockSize);
		coeffsPi.Fill(Pi);
		return Pi;
	}

	SparseMatrix ReduceSparsity(const SparseMatrix& A, const HybridAlgebraicMesh& mesh)
	{
		ThreadLocalCoeffs coeffs;
		#pragma omp parallel for
		for (BigNumber coarseFaceNumber = 0; coarseFaceNumber < mesh.CoarseFaces.size(); ++coarseFaceNumber)
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
						coeffs.Local().Add(it.row(), it.col(), it.value());
				}
			}
		}

		SparseMatrix A2(A.rows(), A.cols());
		coeffs.Fill(A2);

		cout << "A.nonZeros()=" << A.nonZeros() << ", A2.nonZeros()=" << A2.nonZeros() << endl;

		return A2;
	}

public:
	void OnStartSetup() override
	{
		cout << "\t\tk = " << this->PolynomialDegree() << endl;
		cout << "\t\tMesh                : " << this->A_T_T->rows() / _cellBlockSize << " elements, " << this->A_T_F->cols() / _faceBlockSize << " faces";
		if (!this->IsFinestLevel())
		{
			UncondensedLevel* fine = dynamic_cast<UncondensedLevel*>(this->FinerLevel);
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
		if (this->CoarserLevel->ComesFrom == CoarseningType::P)
		{
			// nothing to do
		}
		else
		{
			double scalingFactor = 1.0;
			//scalingFactor = 1.0 / 4.0;
			R = (scalingFactor * P.transpose()).eval();
			//R = scalingFactor * P.transpose();
		}
	}

	Vector Prolong(Vector& vectorOnTheCoarserLevel) override
	{
		if (this->CoarserLevel->ComesFrom == CoarseningType::P)
		{
			// This is possible only if the basis is hierarchical
			UncondensedLevel* coarseLevel = dynamic_cast<UncondensedLevel*>(CoarserLevel);
			auto nHigherDegreeUnknowns = this->_faceBlockSize;
			auto nLowerDegreeUnknowns = coarseLevel->_faceBlockSize;
			auto nFaces = vectorOnTheCoarserLevel.rows() / nLowerDegreeUnknowns;
			Vector vectorOnThisLevel = Vector::Zero(nFaces * nHigherDegreeUnknowns);
			for (BigNumber i = 0; i < nFaces; i++)
				vectorOnThisLevel.segment(i*nHigherDegreeUnknowns, nLowerDegreeUnknowns) = vectorOnTheCoarserLevel.segment(i*nLowerDegreeUnknowns, nLowerDegreeUnknowns);
			return vectorOnThisLevel;
		}
		else
			return Level::Prolong(vectorOnTheCoarserLevel);
	}

	Vector Restrict(Vector& vectorOnThisLevel) override
	{
		if (this->CoarserLevel->ComesFrom == CoarseningType::P)
		{
			// This makes sense only if the basis is hierarchical and orthogonal
			UncondensedLevel* coarseLevel = dynamic_cast<UncondensedLevel*>(CoarserLevel);
			auto nHigherDegreeUnknowns = this->_faceBlockSize;
			auto nLowerDegreeUnknowns = coarseLevel->_faceBlockSize;
			auto nFaces = vectorOnThisLevel.rows() / nHigherDegreeUnknowns;
			Vector vectorOnTheCoarserLevel(nFaces * nLowerDegreeUnknowns);
			for (BigNumber i = 0; i < nFaces; i++)
				vectorOnTheCoarserLevel.segment(i*nLowerDegreeUnknowns, nLowerDegreeUnknowns) = vectorOnThisLevel.segment(i*nHigherDegreeUnknowns, nLowerDegreeUnknowns);
			return vectorOnTheCoarserLevel;
		}
		else
			return Level::Restrict(vectorOnThisLevel);
	}

	Flops ProlongCost() override
	{
		if (this->CoarserLevel->ComesFrom == CoarseningType::P)
			return 0;
		else
			return Level::ProlongCost();
	}

	Flops RestrictCost() override
	{
		if (this->CoarserLevel->ComesFrom == CoarseningType::P)
			return 0;
		else
			return Level::RestrictCost();
	}

	void OnEndSetup() override
	{
		if (this->CoarserLevel)
		{
			UncondensedLevel* coarse = dynamic_cast<UncondensedLevel*>(this->CoarserLevel);
			if (this->CoarserLevel->ComesFrom == CoarseningType::P)
				coarse->SetupOperatorByBlockExtraction();
			else
			{
				coarse->OperatorMatrix = &Ac;
				coarse->A_T_T = &A_T_Tc;
				coarse->A_T_F = &A_T_Fc;
				coarse->A_F_F = ComputeCoarseA_F_F ? &A_F_Fc : nullptr;
			}
		}
	}

	~UncondensedLevel()
	{
		if (!this->IsFinestLevel() && !this->UseGalerkinOperator)
			delete OperatorMatrix;
	}
};