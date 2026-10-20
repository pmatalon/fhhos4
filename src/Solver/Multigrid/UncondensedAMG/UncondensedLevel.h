#pragma once
#include <mutex>
#include <unsupported/Eigen/SparseExtra>
#include "HybridAlgebraicMesh.h"
#include "CoarseningPass.h"
#include "DistributedPass.h"
#include "LocalMatrices.h"
#include "../AggregAMG/AlgebraicMesh.h"
#include "../DownwindOrder.h"
#include "../Level.h"
#include "../../../Utils/SparseMatrixOps.h"
using namespace std;

// The aggregation of the faces of a mesh, as the prolongations need it (UncondensedLevel::BuildProlongation()): whether
// each face is removed (interior to an aggregate), the aggregate of each removed face, and the cells of each face
// (increasing). From HybridAlgebraicMesh, or from the flat coarsening passes (CoarseningPass).
struct FaceAggregation
{
	vector<bool> IsRemoved;
	vector<BigNumber> RemovedFaceAggregate; // -1 for the kept faces
	IndexLists FaceCells;
	BigNumber NCells = 0;
	BigNumber NAggregates = 0;

	FaceAggregation() {}

	// The mesh must be built (its faces' elements) and coarsened
	FaceAggregation(const HybridAlgebraicMesh& mesh)
	{
		BigNumber nFaces = mesh.Faces.size();
		NCells = mesh.Elements.size();
		NAggregates = mesh.CoarseElements.size();
		IsRemoved.resize(nFaces);
		RemovedFaceAggregate.assign(nFaces, -1);
		for (BigNumber f = 0; f < nFaces; ++f)
		{
			IsRemoved[f] = mesh.Faces[f].IsRemovedOnCoarseMesh;
			if (IsRemoved[f])
				RemovedFaceAggregate[f] = mesh.Faces[f].CoarseElements[0]->Number;
		}
		FaceCells.Allocate(nFaces, [&](BigNumber f) { return mesh.Faces[f].Elements.size(); });
		for (BigNumber f = 0; f < nFaces; ++f)
			for (int a = 0; a < FaceCells.Size(f); ++a)
				FaceCells[f][a] = mesh.Faces[f].Elements[a]->Number;
	}
};

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
	// Non-symmetric matrix: the block A_F_T (faces x cells); null if the matrix is symmetric (A_F_T = A_T_F^T)
	const SparseMatrix* A_F_T = nullptr;
private:
	//SparseMatrix Q_T;
	SparseMatrix Q_F;

	SparseMatrix A_T_Tc;
	SparseMatrix A_T_Fc;
	SparseMatrix A_F_Fc;
	SparseMatrix A_F_Tc;
	SparseMatrix _extractedA_F_T;     // A_F_T of a p-coarse level (SetupOperatorByBlockExtraction())
	vector<BigNumber> _downwindOrder; // non-symmetric matrix: the order of the faces in the sweeps of the smoothers
	// On several processes (see FaceMap): the faces of A_T_Fc (the FaceMap of the coarser level), the coarse operator
	// (that of the coarser level), and the operator of a p-coarse level
	shared_ptr<const IndexMap> _coarseFaceMap;
	DistributedMatrix _distributedAc;
	DistributedMatrix _extractedOperator;
	SparseMatrix _extractedA_T_T, _extractedA_T_F;
public:
	SparseMatrix Ac;

	// On several processes (MPI, DistributedOperator set): the local numbers of the faces of the columns of A_T_F, i.e.
	// of the faces of the cells of the process: its own faces (those of the rows of the operator), then the others
	// (ghosts). The cells of A_T_T and A_T_F are all owned. Null on one process.
	shared_ptr<const IndexMap> FaceMap;

	// Set by UncondensedAMG. The same algorithm runs for a symmetric matrix or not; a symmetric one is optimized: the
	// blocks A_T_T, A_T_F, A_F_F and the operator are read from their lower triangular part, the blocks of A_T_T are
	// factorized by Cholesky (LU otherwise), the local Galerkin products (LocalOp) compute the lower triangular parts.
	// Non-symmetric: all the blocks are read in full, A_F_T is given, the coarse operators are P^T A P (restriction P^T).
	bool Symmetric = true;

	// Set by UncondensedAMG: with a non-symmetric matrix, the Gauss-Seidel smoothers follow the downwind order of the
	// operator (DownwindOrder()), in the pre- and post-smoothing
	bool DownwindOrdering = true;

	// Set by UncondensedAMG, to skip the computations that the chosen options don't use:
	// - the blocks A_F_F of the coarse levels, computed at each coarsening step, and the fine one (may be null);
	// - the face prolongation Q_F chained over the coarsening steps.
	// If not computed, the coarse level's A_F_F is null.
	bool ComputeCoarseA_F_F = true;
	bool ComputeQ_F = true;

	// Set by UncondensedAMG: in the face prolongation FaceProlongation, the removed faces that belong to a coarse face
	// take its value instead of the average of the faces of their coarse element (see BuildQ_F()).
	bool ManageAnisotropy = true;

	// Local matrices of the operator, shared by the levels (set by UncondensedAMG). If null, the Galerkin products
	// of the coarsening passes are global sparse products.
	LocalOperator* LocalOp = nullptr;

	// Interpolation of the function 1 on the bases of the cells and faces of this level, one coefficient c per cell and
	// per face (hierarchical bases with a constant first function phi_0: 1 = c phi_0), used to rescale the matrices
	// before the h-coarsening (see CoarsenMesh()). The same on the p-coarse levels (same phi_0). Empty if c = 1
	// everywhere: when it is given so, and on the levels built by h-coarsening, where the constant has coordinate 1 by
	// construction.
	Vector CellInterpOfOne;
	Vector FaceInterpOfOne;

public:
	// The coarsening passes with the default options run on flat arrays (CoarseningPass, see FlatPassesApply()). False
	// runs those of HybridAlgebraicMesh, which give the same results and implement the other options: the reference
	// against which the tests check the flat passes.
	inline static bool UseFlatPasses = true;

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
		if (FaceMap)
			return FaceMap->GlobalSize * _faceBlockSize; // of all the processes
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

		if (this->IsDistributed())
		{
			CoarsenMeshDistributed(coarseningStgy, faceCoarseningStgy, requestedCoarseningFactor, coarsestPossibleMeshReached);
			return;
		}

		if (!this->OperatorMatrix)
		{
			if (this->UseGalerkinOperator)
				ComputeGalerkinOperator();
			else
				SetupDiscretizedOperator();
		}

		if (CellInterpOfOne.rows() == 0)
		{
			CoarsenInPasses(coarseningStgy, faceCoarseningStgy, requestedCoarseningFactor, coarsestPossibleMeshReached);
			return;
		}

		// Correction of a restrictive assumption of the paper (D. A. Di Pietro, F. Hülsemann, P. Matalon, P. Mycek,
		// U. Rüde, "Algebraic multigrid preconditioner for statically condensed systems arising from lowest-order hybrid
		// discretizations", SISC 2023). Its Section 2 takes the DoFs as values ("one scalar value per cell and per
		// face"), so that the coefficients 1 of Q_T, Q_F (3.3) and of the trace Π^f_c (after (3.7)) transfer the
		// constant functions: a coarse cell or face takes the value of the fine ones it aggregates, the trace of a cell
		// value on its faces is the same value. This only holds if the constant function 1 has the same coordinate c in
		// all the cell and face bases (c: the interpolation of 1 on each cell and face basis, 1 = c phi_0). Otherwise (e.g. orthonormal bases: c = sqrt(|T|) on a cell, sqrt(|F|) on a face),
		// the trace of the constant is c_F/c_T, not 1 (it is the geometric trace M_F^-1 M_FT of the h-multigrid that this
		// prolongation mimics), and aggregating fine functions of different coordinates does not give a constant. The
		// coarse spaces then lose the constants (more iterations), or, if c_F < c_T, the prolongation is amplified at
		// each coarsening pass until it overflows.
		// Correction: a change of coordinates on the first (constant) DoF of each cell and face block, z_0 = y_0 / c, in
		// which the constant has coordinate 1 everywhere, as the paper assumes. The coarsening passes run on D A D
		// (D = diag(c) on the first DoFs, 1 elsewhere), and the prolongation is brought back to the coordinates of this
		// level, P = D_F P_z: the operator of this level, its smoothers and residuals stay in the caller's coordinates
		// (the block smoothers are invariant by this scaling, and P^T A P = P_z^T (D_F A D_F) P_z). The coarse levels are
		// in the coordinates z, where the constants have coordinate 1: they need no rescaling.
		const SparseMatrix* levelA     = this->OperatorMatrix;
		const SparseMatrix* levelA_T_T = this->A_T_T;
		const SparseMatrix* levelA_T_F = this->A_T_F;
		const SparseMatrix* levelA_F_F = this->A_F_F;
		const SparseMatrix* levelA_F_T = this->A_F_T;
		SparseMatrix scaledA     = ScaleConstantDoFs(*levelA,     FaceInterpOfOne, _faceBlockSize, FaceInterpOfOne, _faceBlockSize);
		SparseMatrix scaledA_T_T = ScaleConstantDoFs(*levelA_T_T, CellInterpOfOne, _cellBlockSize, CellInterpOfOne, _cellBlockSize);
		SparseMatrix scaledA_T_F = ScaleConstantDoFs(*levelA_T_F, CellInterpOfOne, _cellBlockSize, FaceInterpOfOne, _faceBlockSize);
		SparseMatrix scaledA_F_F = levelA_F_F ? ScaleConstantDoFs(*levelA_F_F, FaceInterpOfOne, _faceBlockSize, FaceInterpOfOne, _faceBlockSize) : SparseMatrix();
		SparseMatrix scaledA_F_T = levelA_F_T ? ScaleConstantDoFs(*levelA_F_T, FaceInterpOfOne, _faceBlockSize, CellInterpOfOne, _cellBlockSize) : SparseMatrix();
		this->OperatorMatrix = &scaledA;
		this->A_T_T = &scaledA_T_T;
		this->A_T_F = &scaledA_T_F;
		this->A_F_F = levelA_F_F ? &scaledA_F_F : nullptr;
		this->A_F_T = levelA_F_T ? &scaledA_F_T : nullptr;

		CoarsenInPasses(coarseningStgy, faceCoarseningStgy, requestedCoarseningFactor, coarsestPossibleMeshReached);

		this->OperatorMatrix = levelA;
		this->A_T_T = levelA_T_T;
		this->A_T_F = levelA_T_F;
		this->A_F_F = levelA_F_F;
		this->A_F_T = levelA_F_T;
		if (LocalOp && LocalOp->Assembled == &scaledA)
			LocalOp->Assembled = nullptr;
		if (coarsestPossibleMeshReached)
			return;

		ScaleConstantRows(this->P, FaceInterpOfOne, _faceBlockSize);
		if (ComputeQ_F)
			ScaleConstantRows(this->Q_F, FaceInterpOfOne, _faceBlockSize);
	}

	// Copy of A whose coefficient (i, j) is multiplied by d_i d_j, d being the interpolation of 1 of the block (the
	// coordinate c of the constant) on the first DoF of each block of rows (resp. columns), and 1 on the other DoFs
	static SparseMatrix ScaleConstantDoFs(const SparseMatrix& A, const Vector& rowInterpOfOne, int rowBlockSize, const Vector& colInterpOfOne, int colBlockSize)
	{
		SparseMatrix scaled = A;
		BigNumber nRows = scaled.rows();
		#pragma omp parallel for
		for (BigNumber i = 0; i < nRows; ++i)
		{
			double d_i = i % rowBlockSize == 0 ? rowInterpOfOne[i / rowBlockSize] : 1;
			for (SparseMatrix::InnerIterator it(scaled, i); it; ++it)
			{
				BigNumber j = it.col();
				double d_j = j % colBlockSize == 0 ? colInterpOfOne[j / colBlockSize] : 1;
				it.valueRef() *= d_i * d_j;
			}
		}
		return scaled;
	}

	// Multiplies the row of the first DoF of each block by the interpolation of 1 of the block (the coordinate c of the
	// constant)
	static void ScaleConstantRows(SparseMatrix& M, const Vector& interpOfOne, int blockSize)
	{
		BigNumber nBlocks = interpOfOne.rows();
		#pragma omp parallel for
		for (BigNumber b = 0; b < nBlocks; ++b)
		{
			for (SparseMatrix::InnerIterator it(M, b * blockSize); it; ++it)
				it.valueRef() *= interpOfOne[b];
		}
	}

	//------------------------------------------------------------------------------------------------------------------//
	//                                        On several processes (MPI)                                                //
	//------------------------------------------------------------------------------------------------------------------//

	// CoarsenMesh() on several processes: the flat coarsening passes (CoarsenInDistributedPasses()), with the default
	// options only, after the rescaling of the constants as on one process
	void CoarsenMeshDistributed(H_CoarsStgy coarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, double requestedCoarseningFactor, bool& coarsestPossibleMeshReached)
	{
		if (!FlatPassesApply(coarseningStgy, faceCoarseningStgy) || _multigridProlong != UAMGProlongation::ChainedCoarseningProlongations || ComputeCoarseA_F_F || ComputeQ_F)
			Utils::FatalError("UncondensedAMG: on several processes, only the default coarsening runs (pairwise aggregation, interface collapsing, coarsening prolongations 3 to 6, chained as multigrid prolongation).");
		if (!Symmetric)
			Utils::FatalError("UncondensedAMG: non-symmetric matrices are not solved on several processes yet.");

		if (CellInterpOfOne.rows() == 0)
		{
			CoarsenInDistributedPasses(coarseningStgy, requestedCoarseningFactor, coarsestPossibleMeshReached);
			return;
		}

		// The rescaling of CoarsenMesh(): the operator's columns of the ghost faces need their interpolations of 1
		const DistributedMatrix* levelA = this->DistributedOperator;
		const SparseMatrix* levelA_T_T = this->A_T_T;
		const SparseMatrix* levelA_T_F = this->A_T_F;
		Vector ownedFaceInterpOfOne = FaceInterpOfOne.head(FaceMap->NOwned);
		Vector ghostFaceInterpOfOne(levelA->ColMap->NGhosts());
		levelA->ColMap->UpdateGhosts(ownedFaceInterpOfOne.data(), ghostFaceInterpOfOne.data(), 1);
		DistributedMatrix scaledA;
		scaledA.RowMap = levelA->RowMap;
		scaledA.ColMap = levelA->ColMap;
		scaledA.RowBlockSize = levelA->RowBlockSize;
		scaledA.ColBlockSize = levelA->ColBlockSize;
		scaledA.Diag = ScaleConstantDoFs(levelA->Diag, ownedFaceInterpOfOne, _faceBlockSize, ownedFaceInterpOfOne, _faceBlockSize);
		scaledA.Offd = ScaleConstantDoFs(levelA->Offd, ownedFaceInterpOfOne, _faceBlockSize, ghostFaceInterpOfOne, _faceBlockSize);
		SparseMatrix scaledA_T_T = ScaleConstantDoFs(*levelA_T_T, CellInterpOfOne, _cellBlockSize, CellInterpOfOne, _cellBlockSize);
		SparseMatrix scaledA_T_F = ScaleConstantDoFs(*levelA_T_F, CellInterpOfOne, _cellBlockSize, FaceInterpOfOne, _faceBlockSize);
		this->DistributedOperator = &scaledA;
		this->OperatorMatrix = &scaledA.Diag;
		this->A_T_T = &scaledA_T_T;
		this->A_T_F = &scaledA_T_F;

		CoarsenInDistributedPasses(coarseningStgy, requestedCoarseningFactor, coarsestPossibleMeshReached);

		this->DistributedOperator = levelA;
		this->OperatorMatrix = &levelA->Diag;
		this->A_T_T = levelA_T_T;
		this->A_T_F = levelA_T_F;
		if (LocalOp->AssembledDistributed == &scaledA)
			LocalOp->AssembledDistributed = nullptr;
		if (coarsestPossibleMeshReached)
			return;

		ScaleConstantRows(this->DistributedP.Diag, ownedFaceInterpOfOne, _faceBlockSize);
		ScaleConstantRows(this->DistributedP.Offd, ownedFaceInterpOfOne, _faceBlockSize);
	}

	// CoarsenInFlatPasses() on several processes: each process aggregates its cells (see DistributedPass), the coarse
	// faces between processes are numbered by their owners, the local matrices of the operator are those of the cells
	// of the process (their faces: those of FaceMap). The prolongations of the passes are chained on all the local faces
	// (the faces between processes are kept: their rows are injections, known on both sides), then the rows of the
	// faces of the process make the distributed prolongation of the level.
	void CoarsenInDistributedPasses(H_CoarsStgy coarseningStgy, double requestedCoarseningFactor, bool& coarsestPossibleMeshReached)
	{
		int fbs = _faceBlockSize;
		const SparseMatrix* passA_T_T = this->A_T_T;
		const SparseMatrix* passA_T_F = this->A_T_F;
		shared_ptr<const IndexMap> passFaceMap = this->FaceMap;
		SparseMatrix A_T_Tc, A_T_Fc;
		SparseMatrix P; // the chained prolongations, rows of all the local faces of the level
		double actualCoarseningFactor = 0;
		int nCoarsenings = 0;
		bool allNeighbours = coarseningStgy == H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours || coarseningStgy == H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours;

		while (!CoarseningCriteriaReached(coarseningStgy, requestedCoarseningFactor, nCoarsenings, actualCoarseningFactor))
		{
			CoarseningPass pass(_cellBlockSize, fbs, Symmetric, _strongCouplingThreshold);
			pass.Build(*passA_T_F, nullptr);
			DistributedPass distributed(pass, *passFaceMap);
			distributed.FindRemoteCells();
			Utils::Log() << (allNeighbours ? "\tElement agglomeration" : "\tElement pairwise aggregation") << " (on each process)" << endl;
			if (pass.Aggregate(allNeighbours))
			{
				coarsestPossibleMeshReached = true;
				return;
			}
			distributed.AddRemoteAggregates();
			Utils::Log() << "\tInterface collapsing" << endl;
			if (pass.CollapseInterfaces())
			{
				coarsestPossibleMeshReached = true;
				return;
			}
			distributed.NumberCoarseFaces();

			// Local matrices of the operator, if they don't come from the previous coarsening pass: of the cells of the
			// process, from the rows of its faces (those of the other processes fetched)
			if (nCoarsenings == 0 && LocalOp->AssembledDistributed != this->DistributedOperator)
			{
				SparseMatrix S = DistributedPass::LocalFaceRows(*this->DistributedOperator, *passFaceMap, fbs);
				LocalOp->Current.Decompose(S, pass.CellFaces, pass.FaceCells, fbs, Symmetric);
			}

			SparseMatrix passP, nextA_T_Tc, nextA_T_Fc;
			pass.ComputeProlongation(_coarseningProlong, _faceProlong, *passA_T_T, *passA_T_F, LocalOp->Current, passP, nextA_T_Tc, nextA_T_Fc);
			LocalOp->Next.GalerkinProduct(LocalOp->Current, pass.AggregateCells, pass.AggregateCoarseFaces, pass.NCoarseFaces, pass.FaceRemoved, pass.FaceCoarseFace, passP);
			LocalOp->Current.Swap(LocalOp->Next);
			LocalOp->Assembled = nullptr;
			LocalOp->AssembledDistributed = nullptr;

			// The coarse faces renumbered as in CoarseFaceMap (owned, then ghosts)
			const vector<BigNumber>& position = distributed.CoarseFacePosition;
			BigNumber nCoarseFaces = distributed.CoarseFaceMap->NLocal();
			auto renumber = [&position](BigNumber cf) { return position[cf]; };
			passP = DistributedMatrix::RenumberColumns(passP, renumber, nCoarseFaces * fbs, fbs);
			nextA_T_Fc = DistributedMatrix::RenumberColumns(nextA_T_Fc, renumber, nCoarseFaces * fbs, fbs);
			for (BigNumber& f : LocalOp->Current.Faces)
				f = position[f];

			if (nCoarsenings == 0)
				P = std::move(passP);
			else
				P = SparseMatrixOps::Multiply(P, passP);

			actualCoarseningFactor = (double)this->FaceMap->GlobalSize / (double)distributed.CoarseFaceMap->GlobalSize;
			Utils::Log() << "\t\tCoarsening factor = " << actualCoarseningFactor << endl;
			if (actualCoarseningFactor >= 10)
				Utils::Warning("The coarsening seems a little too strong...");

			A_T_Tc = std::move(nextA_T_Tc);
			A_T_Fc = std::move(nextA_T_Fc);
			passA_T_T = &A_T_Tc;
			passA_T_F = &A_T_Fc;
			passFaceMap = distributed.CoarseFaceMap;
			nCoarsenings++;
		}

		this->A_T_Tc = std::move(A_T_Tc);
		this->A_T_Fc = std::move(A_T_Fc);
		_coarseFaceMap = passFaceMap;
		// The operator of the level, assembled once (the rows of the coarse faces of other processes added to theirs)
		SparseMatrix localRows = LocalOp->Current.Assemble(_coarseFaceMap->NLocal());
		_distributedAc = DistributedMatrix::AssembleFromLocalRows(_coarseFaceMap, localRows, fbs);
		LocalOp->AssembledDistributed = &_distributedAc;
		// The prolongation: the rows of the faces of the process
		auto rows = make_shared<const IndexMap>(this->FaceMap->WithoutGhosts());
		this->DistributedP = DistributedMatrix::FromLocalColumns(rows, _coarseFaceMap, SparseMatrix(P.topRows(rows->NOwned * fbs)), fbs, fbs);
	}

	// The coarsening passes of the paper, until the requested coarsening factor is reached
	void CoarsenInPasses(H_CoarsStgy coarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, double requestedCoarseningFactor, bool& coarsestPossibleMeshReached)
	{
		if (FlatPassesApply(coarseningStgy, faceCoarseningStgy))
		{
			CoarsenInFlatPasses(coarseningStgy, requestedCoarseningFactor, coarsestPossibleMeshReached);
			return;
		}

		SparseMatrix *auxP, *auxQ_F, *auxSchur;
		HybridAlgebraicMesh *mesh, *coarseMesh;
		const SparseMatrix* schur = this->OperatorMatrix;
		double actualCoarseningFactor = 0;
		double nCoarsenings = 0;

		HybridAlgebraicMesh initialFineMesh(A_T_T, A_T_F, A_F_F, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold);
		initialFineMesh.A_F_T = A_F_T;
		initialFineMesh.Symmetric = Symmetric;
		mesh = &initialFineMesh;

		while (!CoarseningCriteriaReached(coarseningStgy, requestedCoarseningFactor, nCoarsenings, actualCoarseningFactor))
		{
			// Coarsening
			std::tie(coarseMesh, auxP, auxQ_F, auxSchur, coarsestPossibleMeshReached) = Coarsen(*mesh, schur, coarseningStgy, faceCoarseningStgy);
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
			Utils::Log() << "\t\tCoarsening factor = " << actualCoarseningFactor << endl;
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
				if (mesh->A_F_T) delete mesh->A_F_T;
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
			if (coarseMesh->A_F_T)
				this->A_F_Tc = std::move(*coarseMesh->A_F_T);
			if (LocalGalerkinProducts())
			{
				// The operator of the level is only assembled here (for the smoothers and the coarse solver)
				this->Ac = LocalOp->Current.Assemble(this->A_T_Fc.cols() / _faceBlockSize);
				LocalOp->Assembled = &this->Ac;
			}
			else
				this->Ac = std::move(*schur);
		}
		else
		{
			// Multigrid prolongation
			SparseMatrix Q_T = BuildQ_T(initialFineMesh);
			coarseMesh->BuildElementFaces(); // the coarse cells and faces are used by the reconstruction (Theta)
			SparseMatrix* P = BuildProlongation(_multigridProlong, initialFineMesh, this->OperatorMatrix, *coarseMesh, Q_T, &Q_F);
			SetProlongation(P);

			this->A_T_Tc = std::move(*coarseMesh->A_T_T);
			SparseMatrix Pt = SparseMatrixOps::Transpose(this->P);
			this->A_T_Fc = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(SparseMatrixOps::Transpose(Q_T), *this->A_T_F), this->P);
			if (!Symmetric)
				this->A_F_Tc = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, *this->A_F_T), Q_T);
			if (ComputeCoarseA_F_F)
				this->A_F_Fc = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, *this->A_F_F), this->P);
			this->Ac     = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, *this->OperatorMatrix), this->P);
		}
		if (coarseMesh->A_F_T)
			delete coarseMesh->A_F_T;
		delete coarseMesh;
		delete schur;
	}

	// The coarsening passes run on flat arrays (CoarseningPass), with the same results, with the local Galerkin products
	// (interface collapsing, coarsening prolongations 3 to 6), the aggregations of the cells, and the multigrid
	// prolongations other than FindInteriorThatReconstructs (which needs HybridAlgebraicMesh)
	bool FlatPassesApply(H_CoarsStgy coarseningStgy, FaceCoarseningStrategy faceCoarseningStgy) const
	{
		return UseFlatPasses && LocalGalerkinProducts()
			&& _multigridProlong != UAMGProlongation::FindInteriorThatReconstructs && _multigridProlong != UAMGProlongation::ReconstructionTranspose2Steps
			&& (coarseningStgy == H_CoarsStgy::MultiplePairwiseAggregation || coarseningStgy == H_CoarsStgy::DoublePairwiseAggregation
				|| coarseningStgy == H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours || coarseningStgy == H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours)
			&& faceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsing;
	}

	void CoarsenInFlatPasses(H_CoarsStgy coarseningStgy, double requestedCoarseningFactor, bool& coarsestPossibleMeshReached)
	{
		// The blocks of the current pass: those of the level, then those computed by the previous pass
		const SparseMatrix* schur = this->OperatorMatrix;
		const SparseMatrix* passA_T_T = this->A_T_T;
		const SparseMatrix* passA_T_F = this->A_T_F;
		const SparseMatrix* passA_F_T = this->A_F_T;
		const SparseMatrix* passA_F_F = this->A_F_F;
		SparseMatrix A_T_Tc, A_T_Fc, A_F_Tc, A_F_Fc;
		double actualCoarseningFactor = 0;
		int nCoarsenings = 0;

		// Multigrid prolongation other than the chained coarsening prolongations: built on the fine mesh of the level and
		// the aggregation composed over the passes. Composition: the current aggregate of each fine cell; for each fine
		// face, its current face if it is kept (>= 0), -1 - its aggregate if it is removed. The cells of the fine faces.
		bool chained = _multigridProlong == UAMGProlongation::ChainedCoarseningProlongations;
		vector<BigNumber> fineCellAggregate, fineFaceState;
		IndexLists fineFaceCells;

		while (!CoarseningCriteriaReached(coarseningStgy, requestedCoarseningFactor, nCoarsenings, actualCoarseningFactor))
		{
			CoarseningPass pass(_cellBlockSize, _faceBlockSize, Symmetric, _strongCouplingThreshold);
			pass.Build(*passA_T_F, passA_F_T);
			bool allNeighbours = coarseningStgy == H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours || coarseningStgy == H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours;
			Utils::Log() << (allNeighbours ? "\tElement agglomeration" : "\tElement pairwise aggregation") << endl;
			if (pass.Aggregate(allNeighbours))
			{
				coarsestPossibleMeshReached = true;
				return;
			}
			Utils::Log() << "\tInterface collapsing" << endl;
			if (pass.CollapseInterfaces())
			{
				coarsestPossibleMeshReached = true;
				return;
			}

			// Local matrices of the operator, if they don't come from the previous coarsening pass
			if (nCoarsenings == 0 && LocalOp->Assembled != schur)
			{
				LocalOp->Current.Decompose(*schur, pass.CellFaces, pass.FaceCells, _faceBlockSize, Symmetric);
				LocalOp->Assembled = schur;
			}

			SparseMatrix passP, nextA_T_Tc, nextA_T_Fc, nextA_F_Tc, nextA_F_Fc;
			pass.ComputeProlongation(_coarseningProlong, _faceProlong, *passA_T_T, *passA_T_F, LocalOp->Current, passP, nextA_T_Tc, nextA_T_Fc);
			LocalOp->Next.GalerkinProduct(LocalOp->Current, pass.AggregateCells, pass.AggregateCoarseFaces, pass.NCoarseFaces, pass.FaceRemoved, pass.FaceCoarseFace, passP);
			LocalOp->Current.Swap(LocalOp->Next);
			LocalOp->Assembled = nullptr;
			// The coarse blocks A_F_T (non-symmetric matrix) and A_F_F (if used), as in Coarsen()
			if (!Symmetric || ComputeCoarseA_F_F)
			{
				SparseMatrix Pt = SparseMatrixOps::Transpose(passP);
				if (!Symmetric)
					nextA_F_Tc = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, *passA_F_T), pass.Q_T());
				if (ComputeCoarseA_F_F)
					nextA_F_Fc = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, Symmetric ? SparseMatrixOps::FullFromLower(*passA_F_F) : *passA_F_F), passP);
			}

			// Global prolongation operator, and the face prolongation Q_F chained over the passes (if used)
			if (ComputeQ_F)
			{
				if (nCoarsenings == 0)
					this->Q_F = pass.Q_F(_faceProlong);
				else
					this->Q_F = SparseMatrixOps::Multiply(this->Q_F, pass.Q_F(_faceProlong));
			}
			if (chained)
			{
				if (nCoarsenings == 0)
					this->P = std::move(passP);
				else
					this->P = SparseMatrixOps::Multiply(this->P, passP);
			}
			else
			{
				if (nCoarsenings == 0)
				{
					fineCellAggregate = pass.CellAggregate;
					fineFaceCells = pass.FaceCells;
					fineFaceState.resize(pass.NFaces);
					#pragma omp parallel for
					for (BigNumber f = 0; f < pass.NFaces; ++f)
						fineFaceState[f] = f;
				}
				else
				{
					#pragma omp parallel for
					for (BigNumber e = 0; e < (BigNumber)fineCellAggregate.size(); ++e)
						fineCellAggregate[e] = pass.CellAggregate[fineCellAggregate[e]];
				}
				#pragma omp parallel for
				for (BigNumber f = 0; f < (BigNumber)fineFaceState.size(); ++f)
				{
					BigNumber state = fineFaceState[f];
					if (state < 0) // removed in the cell -1 - state of this pass
						fineFaceState[f] = -1 - pass.CellAggregate[-1 - state];
					else if (pass.FaceRemoved[state])
						fineFaceState[f] = -1 - pass.CellAggregate[pass.FaceCells[state][0]];
					else
						fineFaceState[f] = pass.FaceCoarseFace[state];
				}
			}

			actualCoarseningFactor = (double)this->A_T_F->cols() / (double)nextA_T_Fc.cols();
			Utils::Log() << "\t\tCoarsening factor = " << actualCoarseningFactor << endl;
			if (actualCoarseningFactor >= 10)
				Utils::Warning("The coarsening seems a little too strong...");

			A_T_Tc = std::move(nextA_T_Tc);
			A_T_Fc = std::move(nextA_T_Fc);
			A_F_Tc = std::move(nextA_F_Tc);
			A_F_Fc = std::move(nextA_F_Fc);
			passA_T_T = &A_T_Tc;
			passA_T_F = &A_T_Fc;
			passA_F_T = Symmetric ? nullptr : &A_F_Tc;
			passA_F_F = ComputeCoarseA_F_F ? &A_F_Fc : nullptr;
			nCoarsenings++;
		}

		if (!chained)
		{
			// Multigrid prolongation, on the fine mesh of the level (as at the end of CoarsenInPasses())
			FaceAggregation faces;
			faces.NCells = fineCellAggregate.size();
			faces.NAggregates = A_T_Tc.rows() / _cellBlockSize;
			faces.IsRemoved.resize(fineFaceState.size());
			faces.RemovedFaceAggregate.resize(fineFaceState.size());
			for (BigNumber f = 0; f < (BigNumber)fineFaceState.size(); ++f)
			{
				faces.IsRemoved[f] = fineFaceState[f] < 0;
				faces.RemovedFaceAggregate[f] = fineFaceState[f] < 0 ? -1 - fineFaceState[f] : -1;
			}
			faces.FaceCells = std::move(fineFaceCells);
			SparseMatrix Q_T(faces.NCells * _cellBlockSize, faces.NAggregates * _cellBlockSize);
			Q_T.resizeNonZeros(faces.NCells * _cellBlockSize);
			for (BigNumber i = 0; i < faces.NCells * _cellBlockSize; ++i)
			{
				Q_T.outerIndexPtr()[i + 1] = i + 1;
				Q_T.innerIndexPtr()[i] = fineCellAggregate[i / _cellBlockSize] * _cellBlockSize + i % _cellBlockSize;
				Q_T.valuePtr()[i] = 1;
			}
			Q_T.outerIndexPtr()[0] = 0;
			HybridAlgebraicMesh coarseMesh(&A_T_Tc, &A_T_Fc, ComputeCoarseA_F_F ? &A_F_Fc : nullptr, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold);
			coarseMesh.Symmetric = Symmetric;
			coarseMesh.A_F_T = Symmetric ? nullptr : &A_F_Tc;
			coarseMesh.BuildElementFaces(); // the coarse cells and faces are used by the reconstruction (Theta)
			SparseMatrix* P = BuildProlongation(_multigridProlong, faces, *this->A_T_T, *this->A_T_F, this->A_F_F, this->A_F_T, this->OperatorMatrix, coarseMesh, Q_T, &Q_F);
			SetProlongation(P);

			this->A_T_Tc = std::move(A_T_Tc);
			SparseMatrix Pt = SparseMatrixOps::Transpose(this->P);
			this->A_T_Fc = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(SparseMatrixOps::Transpose(Q_T), *this->A_T_F), this->P);
			if (!Symmetric)
				this->A_F_Tc = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, *this->A_F_T), Q_T);
			if (ComputeCoarseA_F_F)
				this->A_F_Fc = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, *this->A_F_F), this->P);
			this->Ac     = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, *this->OperatorMatrix), this->P);
			return;
		}

		this->A_T_Tc = std::move(A_T_Tc);
		this->A_T_Fc = std::move(A_T_Fc);
		if (!Symmetric)
			this->A_F_Tc = std::move(A_F_Tc);
		if (ComputeCoarseA_F_F)
			this->A_F_Fc = std::move(A_F_Fc);
		// The operator of the level is only assembled here (for the smoothers and the coarse solver)
		this->Ac = LocalOp->Current.Assemble(this->A_T_Fc.cols() / _faceBlockSize);
		LocalOp->Assembled = &this->Ac;
	}

	// The multigrid prolongation of the level, from BuildProlongation(), which returns the face prolongation Q_F itself
	// for FaceProlongation: Q_F is then copied (the coarse level may need it), not moved and deleted
	void SetProlongation(SparseMatrix* P)
	{
		if (P == &this->Q_F)
			this->P = this->Q_F;
		else
		{
			this->P = std::move(*P);
			delete P;
		}
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




	// The Galerkin products of the coarsening passes are computed locally (see LocalMatrices) if the coarsening
	// prolongation sends the faces of a cell only to the coarse faces of its aggregate: with the face prolongations
	// with interface collapsing, and the coarsening prolongations made of the face prolongation on the kept faces.
	bool LocalGalerkinProducts() const
	{
		bool interfaceCollapsing = _faceProlong == UAMGFaceProlongation::BoundaryAggregatesInteriorAverage || _faceProlong == UAMGFaceProlongation::BoundaryAggregatesInteriorZero;
		bool localProlongation = _coarseningProlong == UAMGProlongation::FaceProlongation || _coarseningProlong == UAMGProlongation::FaceProlongationAndInteriorSmoothing
			|| _coarseningProlong == UAMGProlongation::ReconstructTraceOrInject || _coarseningProlong == UAMGProlongation::ReconstructSmoothedTraceOrInject;
		return LocalOp && interfaceCollapsing && localProlongation;
	}

	// Returns <coarseMesh, P, Q_F, schurc, coarsestPossibleMeshReached>.
	// With the local Galerkin products, schur is null after the first coarsening pass of the level, and schurc is
	// not assembled (null): the local matrices of the operators are in LocalOp.
	tuple<HybridAlgebraicMesh*, SparseMatrix*, SparseMatrix*, SparseMatrix*, bool> Coarsen(HybridAlgebraicMesh& mesh, const SparseMatrix* schur, H_CoarsStgy elemCoarseningStgy, FaceCoarseningStrategy faceCoarseningStgy)
	{
		bool local = LocalGalerkinProducts();

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

		// Local matrices of the operator, if they don't come from the previous coarsening pass
		if (local && schur && LocalOp->Assembled != schur)
		{
			LocalOp->Current.Decompose(*schur, mesh, _faceBlockSize, Symmetric);
			LocalOp->Assembled = schur;
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
			skeleton.Build(*schur);
			skeleton.PairWiseAggregate(coarsestPossibleMeshReached);
			if (coarsestPossibleMeshReached)
				return { nullptr, nullptr, nullptr, nullptr, coarsestPossibleMeshReached };
			Q_F = new SparseMatrix(BuildQ_F_AllAggregated(skeleton));
		}
		else
			Utils::FatalError("Unmanaged -face-prolong");

		// Intermediate coarse operators
		// The products are evaluated from left to right, and the symmetric matrices are given by their lower triangular part
		// (non-symmetric: in full)
		auto full = [this](const SparseMatrix& M) { return Symmetric ? SparseMatrixOps::FullFromLower(M) : M; };
		SparseMatrix Q_Tt = SparseMatrixOps::Transpose(Q_T);
		SparseMatrix Q_Tt_A_T_F = SparseMatrixOps::Multiply(Q_Tt, *mesh.A_T_F);
		SparseMatrix* A_T_Tc = nullptr;
		SparseMatrix* A_T_Fc_tmp = nullptr;
		if (!onlyFacesUsed)
		{
			A_T_Tc     = new SparseMatrix(SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Q_Tt, full(*mesh.A_T_T)), Q_T));
			A_T_Fc_tmp = new SparseMatrix(SparseMatrixOps::Multiply(Q_Tt_A_T_F, *Q_F));
		}

		// The prolongation only uses its matrices and the faces of its elements (reconstruction, Theta())
		HybridAlgebraicMesh auxCoarseMesh(A_T_Tc, A_T_Fc_tmp, nullptr, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold);
		auxCoarseMesh.Symmetric = Symmetric;
		auxCoarseMesh.BuildElementFaces();

		/*ExportMatrix(Q_T, "Q_T", 0);
		ExportMatrix(*Q_F, "Q_F", 0);
		ExportMatrix(*A_T_Tc, "A_T_Tc", 0);
		ExportMatrix(*A_T_Fc_tmp, "A_T_Fc", 0);*/

		// Multigrid prolongation
		SparseMatrix* P = BuildProlongation(this->_coarseningProlong, mesh, schur, auxCoarseMesh, Q_T, Q_F, local ? &LocalOp->Current : nullptr);

		SparseMatrix* A_T_Fc = new SparseMatrix(SparseMatrixOps::Multiply(Q_Tt_A_T_F, *P)); // Kills -prolong 1 or 2 because P is then very dense
		//SparseMatrix* A_T_Fc = A_T_Fc_tmp;

		SparseMatrix Pt; // for A_F_Fc, the global Galerkin product, and A_F_Tc (non-symmetric matrix)
		if (ComputeCoarseA_F_F || !local || !Symmetric)
			Pt = SparseMatrixOps::Transpose(*P);
		SparseMatrix* A_F_Fc = nullptr;
		if (ComputeCoarseA_F_F)
			A_F_Fc = new SparseMatrix(SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, full(*mesh.A_F_F)), *P));
		SparseMatrix* schurc = nullptr;
		if (local)
		{
			LocalOp->Next.GalerkinProduct(LocalOp->Current, mesh, *P);
			LocalOp->Current.Swap(LocalOp->Next);
			LocalOp->Assembled = nullptr;
		}
		else
			schurc = new SparseMatrix(SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, full(*schur)), *P));

		HybridAlgebraicMesh* coarseMesh = new HybridAlgebraicMesh(A_T_Tc, A_T_Fc, A_F_Fc, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold);
		coarseMesh->Symmetric = Symmetric;
		if (!Symmetric)
			coarseMesh->A_F_T = new SparseMatrix(SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, *mesh.A_F_T), Q_T));
		
		return { coarseMesh, P, Q_F, schurc, coarsestPossibleMeshReached };
	}



	// The block Jacobi smoothing of the prolongations 4 and 6 uses the rows of the removed faces of the operator, schur,
	// or, if they are given, the local matrices of the operator (schur is then not used).
	SparseMatrix* BuildProlongation(UAMGProlongation prolong, HybridAlgebraicMesh& mesh, const SparseMatrix* schur, HybridAlgebraicMesh& coarseMesh,
									const SparseMatrix& Q_T, SparseMatrix* Q_F, const LocalMatrices* localSchur = nullptr)
	{
		if (prolong == UAMGProlongation::FindInteriorThatReconstructs)
			return BuildProlongationFindInteriorThatReconstructs(mesh, coarseMesh, Q_T, Q_F);
		return BuildProlongation(prolong, FaceAggregation(mesh), *mesh.A_T_T, *mesh.A_T_F, mesh.A_F_F, mesh.A_F_T, schur, coarseMesh, Q_T, Q_F, localSchur);
	}

	// Same, the mesh given by its blocks (A_F_F, A_F_T: null if not given) and the aggregation of its faces, for the
	// prolongations that only need these (all but FindInteriorThatReconstructs)
	SparseMatrix* BuildProlongation(UAMGProlongation prolong, const FaceAggregation& faces, const SparseMatrix& meshA_T_T, const SparseMatrix& meshA_T_F,
	                                const SparseMatrix* meshA_F_F, const SparseMatrix* meshA_F_T, const SparseMatrix* schur, HybridAlgebraicMesh& coarseMesh,
	                                const SparseMatrix& Q_T, SparseMatrix* Q_F, const LocalMatrices* localSchur = nullptr)
	{
		SparseMatrix* P = nullptr;
		if (prolong == UAMGProlongation::ReconstructionTrace) // 1
		{
			// Theta: reconstruction from the coarse faces to the coarse cells
			SparseMatrix Theta = coarseMesh.Theta();
			// Pi: average on both sides of each face
			SparseMatrix Pi = BuildTrace(faces);

			P = new SparseMatrix(SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pi, Q_T), Theta));
		}
		else if (prolong == UAMGProlongation::FaceProlongation) // 3
		{
			// -g 1 -prolong 3 -face-prolong 3 -cs z
			P = Q_F;
		}
		else if (prolong == UAMGProlongation::FaceProlongationAndInteriorSmoothing) // 4
		{
			const vector<bool>& isRemoved = faces.IsRemoved;

			// Smoothing (only the rows of the removed faces are kept)
			BlockJacobi blockJacobi(_faceBlockSize, 2.0 / 3.0);
			SparseMatrix removedRows;
			SetupBlockJacobiOnRemovedFaces(blockJacobi, isRemoved, schur, localSchur, removedRows);
			SparseMatrix J = blockJacobi.IterationMatrix(isRemoved);

			SparseMatrix smoothedQ_F = SparseMatrixOps::Multiply(J, *Q_F);

			P = new SparseMatrix(RemovedOrKeptFaceRows(isRemoved, smoothedQ_F, *Q_F));
		}
		else if (prolong == UAMGProlongation::ReconstructTraceOrInject) // 5
		{
			SparseMatrix Theta = coarseMesh.Theta();   // Reconstruct
			SparseMatrix Pi = BuildCoarseTraceOnFineRemovedFaces(faces); // Trace

			SparseMatrix ReconstructAndTrace = SparseMatrixOps::Multiply(Pi, Theta);

			P = new SparseMatrix(RemovedOrKeptFaceRows(faces.IsRemoved, ReconstructAndTrace, *Q_F));
		}
		else if (prolong == UAMGProlongation::ReconstructSmoothedTraceOrInject) // 6
		{
			const vector<bool>& isRemoved = faces.IsRemoved;

			SparseMatrix Theta = coarseMesh.Theta();   // Reconstruct
			SparseMatrix Pi = BuildCoarseTraceOnFineRemovedFaces(faces); // Trace

			SparseMatrix ReconstructAndTrace = SparseMatrixOps::Multiply(Pi, Theta);

			SparseMatrix ReconstructTraceOrInject = RemovedOrKeptFaceRows(isRemoved, ReconstructAndTrace, *Q_F);

			// Smoothing (only the rows of the removed faces are kept: the matrix J is not assembled for the others)
			BlockJacobi blockJacobi(_faceBlockSize, 2.0/3.0);
			SparseMatrix removedRows;
			SetupBlockJacobiOnRemovedFaces(blockJacobi, isRemoved, schur, localSchur, removedRows);
			SparseMatrix J = blockJacobi.IterationMatrix(isRemoved);

			SparseMatrix ReconstructAndSmoothedTrace = SparseMatrixOps::Multiply(J, ReconstructTraceOrInject);

			P = new SparseMatrix(RemovedOrKeptFaceRows(isRemoved, ReconstructAndSmoothedTrace, *Q_F));
		}
		else if (prolong == UAMGProlongation::HighOrder) // 8
		{
			SparseMatrix Theta = coarseMesh.Theta();
			SparseMatrix Pi = BuildHighOrderTraceOnRemovedFaces(faces, meshA_T_T, meshA_T_F, meshA_F_F, meshA_F_T);

			SparseMatrix ReconstructAndTrace1 = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pi, Q_T), Theta);

			P = new SparseMatrix(RemovedOrKeptFaceRows(faces.IsRemoved, ReconstructAndTrace1, *Q_F));
		}
		else if (prolong == UAMGProlongation::ReconstructionTranspose2Steps) // 9
		{
			Utils::FatalError("The prolongation ReconstructionTranspose2Steps is not implemented.");
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

	// Prolongation FindInteriorThatReconstructs (7): on the removed faces, the values whose reconstruction in the fine
	// cells best fits (least squares) the reconstruction in their coarse cell
	SparseMatrix* BuildProlongationFindInteriorThatReconstructs(HybridAlgebraicMesh& mesh, HybridAlgebraicMesh& coarseMesh, const SparseMatrix& Q_T, SparseMatrix* Q_F)
	{
		SparseMatrix* P;
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
			DenseMatrix Theta_Tc = Symmetric ? DenseMatrix(-A_Tc_Tc.llt().solve(A_Tc_F)) : DenseMatrix(-A_Tc_Tc.partialPivLu().solve(A_Tc_F));

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
				DenseMatrix Theta_Tf = Symmetric ? DenseMatrix(-A_Tf_Tf.llt().solve(A_Tf_F)) : DenseMatrix(-A_Tf_Tf.partialPivLu().solve(A_Tf_F));
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
		return P;
	}

	// Block Jacobi for the rows of the removed faces only: from the operator schur, or from its local matrices, whose
	// rows of the removed faces are assembled in removedRows (which must live as long as blockJacobi is used)
	void SetupBlockJacobiOnRemovedFaces(BlockJacobi& blockJacobi, const vector<bool>& isRemoved, const SparseMatrix* schur, const LocalMatrices* localSchur, SparseMatrix& removedRows)
	{
		if (localSchur)
		{
			removedRows = localSchur->Assemble(isRemoved.size(), &isRemoved);
			blockJacobi.SetupForIterationMatrix(removedRows, isRemoved);
		}
		else
			blockJacobi.SetupForIterationMatrix(*schur, isRemoved);
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
		if (fine->IsDistributed())
		{
			SetupOperatorByBlockExtractionDistributed(*fine);
			return;
		}

		auto nElems = fine->A_T_F->rows() / fine->_cellBlockSize;
		auto nFaces = fine->A_T_F->cols() / fine->_faceBlockSize;

		this->OperatorMatrix = ExtractCoarseMatrix(*this->FinerLevel->OperatorMatrix, nFaces, nFaces, fine->_faceBlockSize, fine->_faceBlockSize, this->_faceBlockSize, this->_faceBlockSize);
		this->A_T_T = ExtractCoarseMatrix(*fine->A_T_T, nElems, nElems, fine->_cellBlockSize, fine->_cellBlockSize, this->_cellBlockSize, this->_cellBlockSize);
		this->A_T_F = ExtractCoarseMatrix(*fine->A_T_F, nElems, nFaces, fine->_cellBlockSize, fine->_faceBlockSize, this->_cellBlockSize, this->_faceBlockSize);
		if (fine->A_F_T)
		{
			SparseMatrix* A_F_T = ExtractCoarseMatrix(*fine->A_F_T, nFaces, nElems, fine->_faceBlockSize, fine->_cellBlockSize, this->_faceBlockSize, this->_cellBlockSize);
			_extractedA_F_T = std::move(*A_F_T);
			delete A_F_T;
			this->A_F_T = &_extractedA_F_T;
		}
		// A_F_F: only for the options that use it (see UncondensedAMG::A_F_FNeeded())
		this->A_F_F = nullptr;
		if (ComputeCoarseA_F_F)
		{
			if (!fine->A_F_F)
				Utils::FatalError("UncondensedAMG: the p-coarsening requires the block A_F_F, which has not been computed on the finer level (see UncondensedAMG::A_F_FNeeded()).");
			this->A_F_F = ExtractCoarseMatrix(*fine->A_F_F, nFaces, nFaces, fine->_faceBlockSize, fine->_faceBlockSize, this->_faceBlockSize, this->_faceBlockSize);
		}
	}

	// On several processes: the same blocks of the rows of the process, the same faces (FaceMap, and the columns of the
	// operator)
	void SetupOperatorByBlockExtractionDistributed(const UncondensedLevel& fine)
	{
		if (ComputeCoarseA_F_F || fine.A_F_T)
			Utils::FatalError("UncondensedAMG: on several processes, the p-coarsening runs without the blocks A_F_F and A_F_T.");
		const DistributedMatrix& fineA = *fine.DistributedOperator;
		int ffbs = fine._faceBlockSize, fcbs = fine._cellBlockSize;
		BigNumber nElems = fine.A_T_F->rows() / fcbs;
		BigNumber nOwnedFaces = fineA.RowMap->NOwned, nGhostFaces = fineA.ColMap->NGhosts();
		auto extract = [](SparseMatrix* M) { SparseMatrix copy = std::move(*M); delete M; return copy; };
		_extractedOperator.RowMap = fineA.RowMap;
		_extractedOperator.ColMap = fineA.ColMap;
		_extractedOperator.RowBlockSize = _faceBlockSize;
		_extractedOperator.ColBlockSize = _faceBlockSize;
		_extractedOperator.Diag = extract(ExtractCoarseMatrix(fineA.Diag, nOwnedFaces, nOwnedFaces, ffbs, ffbs, _faceBlockSize, _faceBlockSize));
		_extractedOperator.Offd = extract(ExtractCoarseMatrix(fineA.Offd, nOwnedFaces, nGhostFaces, ffbs, ffbs, _faceBlockSize, _faceBlockSize));
		this->DistributedOperator = &_extractedOperator;
		this->OperatorMatrix = &_extractedOperator.Diag;
		_extractedA_T_T = extract(ExtractCoarseMatrix(*fine.A_T_T, nElems, nElems, fcbs, fcbs, _cellBlockSize, _cellBlockSize));
		_extractedA_T_F = extract(ExtractCoarseMatrix(*fine.A_T_F, nElems, fine.FaceMap->NLocal(), fcbs, ffbs, _cellBlockSize, _faceBlockSize));
		this->A_T_T = &_extractedA_T_T;
		this->A_T_F = &_extractedA_T_F;
		this->A_F_F = nullptr;
		this->A_F_T = nullptr;
		this->FaceMap = fine.FaceMap;
	}

private:
	// The first coarseBlockRows x coarseBlockCols coefficients of each block (hierarchical bases), without those of
	// absolute value <= NonZeroCoefficients::ZeroThreshold. The rows keep the order of their columns.
	static SparseMatrix* ExtractCoarseMatrix(const SparseMatrix& fineMatrix, BigNumber nBlockRows, BigNumber nBlockCols, int fineBlockRows, int fineBlockCols, int coarseBlockRows, int coarseBlockCols)
	{
		return new SparseMatrix(SparseMatrixOps::BuildByRows(nBlockRows * coarseBlockRows, nBlockCols * coarseBlockCols, coarseBlockRows,
			[&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			for (BigNumber r = firstRow; r < lastRow; ++r)
			{
				BigNumber i = r / coarseBlockRows;
				int k = r % coarseBlockRows;
				for (RowMajorSparseMatrix::InnerIterator it(fineMatrix, i*fineBlockRows + k); it; ++it)
				{
					auto j = it.col() / fineBlockCols;
					int l = it.col() - j * fineBlockCols;
					if (l < coarseBlockCols && abs(it.value()) > NonZeroCoefficients::ZeroThreshold)
						chunk.Add(j*coarseBlockCols + l, it.value());
				}
				chunk.EndRow();
			}
		}));
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

	// Cell prolongation Q_T with only one 1 coefficient per row: (3.3a) of the paper. Like Q_F and the traces below, it
	// preserves the constants only if they have the same coordinate in all the bases, a restrictive assumption of the
	// paper that the rescaling of CoarsenMesh() makes hold.
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

	// Face prolongation FaceProlongation: (3.3b-c) of the paper (same assumption as BuildQ_T())
	SparseMatrix BuildQ_F(const HybridAlgebraicMesh& mesh)
	{
		bool enableAnisotropyManagement = false;
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);

		ThreadLocalCoeffs coeffs1;
		#pragma omp parallel for
		for (BigNumber faceNumber = 0; faceNumber < mesh.Faces.size(); ++faceNumber)
		{
			const HybridAlgebraicFace* face = &mesh.Faces[faceNumber];
			if (face->IsRemovedOnCoarseMesh && (!ManageAnisotropy || !face->CoarseFace))
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

	SparseMatrix BuildTrace(const FaceAggregation& faces)
	{
		// Pi: average on both sides of each face. Trace of the constant 1: see BuildCoarseTraceOnFineRemovedFaces().
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;

		BigNumber nFaces = faces.IsRemoved.size();
		ThreadLocalCoeffs coeffsPi;
		#pragma omp parallel for
		for (BigNumber faceNumber = 0; faceNumber < nFaces; ++faceNumber)
		{
			const BigNumber* cells = faces.FaceCells[faceNumber];
			int nCells = faces.FaceCells.Size(faceNumber);
			if (faces.IsRemoved[faceNumber])
				coeffsPi.Local().Add(faceNumber*_faceBlockSize, cells[0]*_cellBlockSize, traceOfConstant);
			else
			{
				assert(nCells > 0);
				for (int q = 0; q < nCells; ++q)
					coeffsPi.Local().Add(faceNumber*_faceBlockSize, cells[q]*_cellBlockSize, 1.0 / nCells*traceOfConstant);
			}
		}
		SparseMatrix Pi = SparseMatrix(nFaces*_faceBlockSize, faces.NCells*_cellBlockSize);
		coeffsPi.Fill(Pi);
		return Pi;
	}

	// Trace Π^f_c of the paper (after (3.7)), from the coarse cells to the fine faces they contain
	SparseMatrix BuildCoarseTraceOnFineRemovedFaces(const FaceAggregation& faces)
	{
		// Restrictive assumption of the paper: the trace of the constant is 1, i.e. the constant has the same coordinate
		// c in the cell and face bases. In general, it is c_F/c_T (geometric trace M_F^-1 M_FT): CoarsenMesh() rescales
		// the matrices so that c = 1 everywhere. Only the constant mode is transferred (k >= 1 with -hp-cs h: the higher
		// modes of the faces are then only set by the smoothing of the prolongation).
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;

		BigNumber nFaces = faces.IsRemoved.size();
		ThreadLocalCoeffs coeffsPi;
		#pragma omp parallel for
		for (BigNumber faceNumber = 0; faceNumber < nFaces; ++faceNumber)
		{
			if (faces.IsRemoved[faceNumber])
				coeffsPi.Local().Add(faceNumber*_faceBlockSize, faces.RemovedFaceAggregate[faceNumber]*_cellBlockSize, traceOfConstant);
		}
		SparseMatrix Pi = SparseMatrix(nFaces*_faceBlockSize, faces.NAggregates*_cellBlockSize);
		coeffsPi.Fill(Pi);
		return Pi;
	}

	// Same assumption as BuildCoarseTraceOnFineRemovedFaces()
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

	// The mesh given by its blocks (A_F_T: null if the matrix is symmetric) and the aggregation of its faces
	SparseMatrix BuildHighOrderTraceOnRemovedFaces(const FaceAggregation& faces, const SparseMatrix& meshA_T_T, const SparseMatrix& meshA_T_F, const SparseMatrix* meshA_F_F, const SparseMatrix* meshA_F_T)
	{
		BigNumber nFaces = faces.IsRemoved.size();
		ThreadLocalCoeffs coeffsPi;
		#pragma omp parallel for
		for (BigNumber faceNumber = 0; faceNumber < nFaces; ++faceNumber)
		{
			if (faces.IsRemoved[faceNumber])
			{
				BigNumber elemNumber = faces.FaceCells[faceNumber][0];
				DenseMatrix faceMass     = meshA_F_F->block(faceNumber * _faceBlockSize, faceNumber * _faceBlockSize, _faceBlockSize, _faceBlockSize);
				DenseMatrix cellFaceMass = meshA_T_F.block(elemNumber * _cellBlockSize, faceNumber * _faceBlockSize, _cellBlockSize, _faceBlockSize);
				// Non-symmetric matrix: the block A_F_T(face, elem) instead of A_T_F(elem, face)^T
				DenseMatrix faceCellMass = meshA_F_T ? DenseMatrix(meshA_F_T->block(faceNumber * _faceBlockSize, elemNumber * _cellBlockSize, _faceBlockSize, _cellBlockSize)) : DenseMatrix(cellFaceMass.transpose());
				DenseMatrix trace = Symmetric ? DenseMatrix(-faceMass.llt().solve(faceCellMass)) : DenseMatrix(-faceMass.partialPivLu().solve(faceCellMass));
				//DenseMatrix trace = - cellFaceMass.transpose();
				coeffsPi.Local().Add(faceNumber*_faceBlockSize, elemNumber*_cellBlockSize, trace);
			}
		}
		SparseMatrix Pi = SparseMatrix(nFaces*_faceBlockSize, faces.NCells*_cellBlockSize);
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

		Utils::Log() << "A.nonZeros()=" << A.nonZeros() << ", A2.nonZeros()=" << A2.nonZeros() << endl;

		return A2;
	}

public:
	void OnStartSetup() override
	{
		Utils::Log() << "\t\tk = " << this->PolynomialDegree() << endl;
		// On several processes: the cells and faces of all of them
		BigNumber nElements = Communicator::Current().Sum((BigNumber)this->A_T_T->rows() / _cellBlockSize);
		BigNumber nFaces = this->NUnknowns() / _faceBlockSize;
		Utils::Log() << "\t\tMesh                : " << nElements << " elements, " << nFaces << " faces";
		if (!this->IsFinestLevel())
		{
			UncondensedLevel* fine = dynamic_cast<UncondensedLevel*>(this->FinerLevel);
			double nFine = fine->NUnknowns();
			double nCoarse = this->NUnknowns();
			Utils::Log() << ", coarsening factor = " << (nFine/nCoarse);
		}
		Utils::Log() << endl;
		SetupDownwindOrder();
	}

	// Non-symmetric matrix: the downwind order of the faces, followed by all the sweeps of the smoothers (pre- and
	// post-smoothing) if some couplings are one-sided (convection), the natural order otherwise
	void SetupDownwindOrder()
	{
		if (Symmetric || !DownwindOrdering || this->IsCoarsestLevel())
			return;
		BigNumber nOneSided, nCycleBreaks;
		_downwindOrder = DownwindOrder(*this->OperatorMatrix, _faceBlockSize, FaceInterpOfOne, _strongCouplingThreshold, nOneSided, nCycleBreaks);
		Utils::Log() << "\t\tDownwind order      : ";
		if (_downwindOrder.empty())
			Utils::Log() << "none (no one-sided coupling): natural order" << endl;
		else
		{
			Utils::Log() << nOneSided << " one-sided couplings, " << nCycleBreaks << " cycles broken" << endl;
			this->PreSmoother->SetSweepOrder(&_downwindOrder);
			this->PostSmoother->SetSweepOrder(&_downwindOrder);
		}
	}

	void SetupProlongation() override
	{}

	void SetupRestriction() override
	{
		if (this->CoarserLevel->ComesFrom == CoarseningType::P)
		{
			// nothing to do
		}
		else if (this->IsDistributed())
			this->DistributedP.SetupTransposed();
		else
			R = SparseMatrixOps::Transpose(P);
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
			{
				coarse->SetupOperatorByBlockExtraction();
				// The lower-degree bases start with the same constant function
				coarse->CellInterpOfOne = this->CellInterpOfOne;
				coarse->FaceInterpOfOne = this->FaceInterpOfOne;
			}
			else if (this->IsDistributed())
			{
				coarse->DistributedOperator = &_distributedAc;
				coarse->OperatorMatrix = &_distributedAc.Diag;
				coarse->A_T_T = &A_T_Tc;
				coarse->A_T_F = &A_T_Fc;
				coarse->A_F_F = nullptr;
				coarse->A_F_T = nullptr;
				coarse->FaceMap = _coarseFaceMap;
			}
			else
			{
				coarse->OperatorMatrix = &Ac;
				coarse->A_T_T = &A_T_Tc;
				coarse->A_T_F = &A_T_Fc;
				coarse->A_F_F = ComputeCoarseA_F_F ? &A_F_Fc : nullptr;
				coarse->A_F_T = Symmetric ? nullptr : &A_F_Tc;
			}
		}
	}

	~UncondensedLevel()
	{
		if (!this->IsFinestLevel() && !this->UseGalerkinOperator)
			delete OperatorMatrix;
	}
};