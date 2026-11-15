// fhhos4's preconditioner for PETSc, PCFHHOS4 (-pc_type fhhos4): the U-AMG of the library for the uncondensed system
// of a hybrid discretization, given as a MATNEST [[A_TT, A_TF], [A_FT, A_FF]] of AIJ blocks (the cells, then the
// faces). One application: one cycle (Krylov::None, MaxIterations 1), i.e. the condensation of the residual, one cycle
// on the condensed system, and the cell unknowns. The default K-cycle is non-linear: PETSc's Krylov method must be
// flexible (fcg, fgmres, gcr). Written in PETSc's style, on PETSc's interface for the implementations of
// preconditioners and the public header of the library, with its inference of the discretization (header-only,
// compiled into this library), so as to move into PETSc (src/ksp/pc/impls) later.
#include <petsc/private/pcimpl.h>
#include <fhhos4>
#include "fhhos4_petsc.h"
#include "Solver/DiscretizationInference.h"
#include <cmath>
#include <new>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

static_assert(std::is_same<PetscScalar, double>::value, "fhhos4: a PETSc with real double-precision scalars");

namespace
{
	using Map = Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor, PetscInt>>;
	// A block for the setups: the map of PETSc's arrays if its index type is the library's, a copy otherwise (by Eigen)
	using LibraryMatrix = std::conditional_t<std::is_same<PetscInt, fhhos4::SparseIndex>::value, Map, fhhos4::SparseMatrix>;

	struct PC_FHHOS4
	{
		fhhos4::Solver Solver;                 // its parameters: the options -pc_fhhos4_*
		PetscBool SymmetricSet = PETSC_FALSE;  // -pc_fhhos4_symmetric given; otherwise MatIsSymmetricKnown() of the operator
		PetscBool Symmetric = PETSC_FALSE;
		PetscBool SetUpSymmetric = PETSC_FALSE; // that of the last setup
		IS Cells = nullptr, Faces = nullptr;    // the rows of the blocks in the operator (MatNestGetISs())
		fhhos4::DistributedBlocks Merged;       // several processes: the blocks merged (MergeOffDiag())
		fhhos4::SparseMatrix A_FT;              // one process, non-symmetric: a copy, read by the solves
	};

	// Whether ok on all the processes
	PetscErrorCode AllOk(MPI_Comm comm, bool ok, bool* all)
	{
		int local = ok, global = 1;
		PetscFunctionBegin;
		PetscCallMPI(MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_LAND, comm));
		*all = global;
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	// An error on all the processes if one of them fails the check (the others would wait in the next collective call)
	PetscErrorCode CheckOnAll(MPI_Comm comm, bool ok, const std::string& message)
	{
		bool all;
		PetscFunctionBegin;
		PetscCall(AllOk(comm, ok, &all));
		PetscCheck(all, comm, PETSC_ERR_ARG_WRONG, "fhhos4: %s", ok ? "an input is wrong on another process (see its error)" : message.c_str());
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	// The block (i, j) of the MATNEST, an AIJ matrix
	PetscErrorCode SubMatrix(Mat nest, PetscInt i, PetscInt j, const char* name, Mat* A)
	{
		PetscBool aij;
		PetscFunctionBegin;
		PetscCall(MatNestGetSubMat(nest, i, j, A));
		PetscCheck(*A, PetscObjectComm((PetscObject)nest), PETSC_ERR_ARG_WRONG, "fhhos4: the block %s of the MATNEST is missing", name);
		PetscCall(PetscObjectBaseTypeCompareAny((PetscObject)*A, &aij, MATSEQAIJ, MATMPIAIJ, ""));
		PetscCheck(aij, PetscObjectComm((PetscObject)*A), PETSC_ERR_SUP, "fhhos4: the block %s of the MATNEST must be an AIJ matrix, not %s (MatConvert())", name, ((PetscObject)*A)->type_name);
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	// The CSR arrays of a SEQAIJ matrix, read in place until Restore()
	struct SeqAIJ
	{
		Mat M = nullptr;
		PetscInt Rows = 0, Cols = 0;
		const PetscInt *RowPointers = nullptr, *Columns = nullptr;
		const PetscScalar* Values = nullptr;

		PetscErrorCode Get(Mat m)
		{
			PetscBool done;
			PetscFunctionBegin;
			M = m;
			PetscCall(MatGetSize(M, &Rows, &Cols));
			PetscCall(MatGetRowIJ(M, 0, PETSC_FALSE, PETSC_FALSE, &Rows, &RowPointers, &Columns, &done));
			PetscCheck(done, PETSC_COMM_SELF, PETSC_ERR_SUP, "fhhos4: MatGetRowIJ() failed");
			PetscCall(MatSeqAIJGetArrayRead(M, &Values));
			PetscFunctionReturn(PETSC_SUCCESS);
		}

		PetscErrorCode Restore()
		{
			PetscBool done;
			PetscFunctionBegin;
			if (M)
			{
				PetscCall(MatRestoreRowIJ(M, 0, PETSC_FALSE, PETSC_FALSE, &Rows, &RowPointers, &Columns, &done));
				PetscCall(MatSeqAIJRestoreArrayRead(M, &Values));
				M = nullptr;
			}
			PetscFunctionReturn(PETSC_SUCCESS);
		}

		Map Eigen() const
		{
			return Map(Rows, Cols, RowPointers ? RowPointers[Rows] : 0, RowPointers, Columns, Values);
		}
	};

	// A block of the operator, an AIJ matrix: on several processes (MATMPIAIJ), in two parts, Diag, Offd, and the global
	// numbers of the columns of Offd (OffdColumns); on one process, Diag alone (MATSEQAIJ, or MATMPIAIJ, whose Diag is
	// then the whole block)
	struct Block
	{
		SeqAIJ Diag, Offd;
		const PetscInt* OffdColumns = nullptr;

		PetscErrorCode Get(Mat A)
		{
			PetscBool mpi;
			PetscFunctionBegin;
			PetscCall(PetscObjectBaseTypeCompare((PetscObject)A, MATMPIAIJ, &mpi));
			if (mpi)
			{
				Mat Ad, Ao;
				PetscCall(MatMPIAIJGetSeqAIJ(A, &Ad, &Ao, &OffdColumns));
				PetscCall(Diag.Get(Ad));
				PetscCall(Offd.Get(Ao));
			}
			else
				PetscCall(Diag.Get(A));
			PetscFunctionReturn(PETSC_SUCCESS);
		}

		PetscErrorCode Restore()
		{
			PetscFunctionBegin;
			PetscCall(Diag.Restore());
			PetscCall(Offd.Restore());
			PetscFunctionReturn(PETSC_SUCCESS);
		}
	};

	// The near null space of the operator on the DoFs of the cells and of the faces of the process, one column per
	// vector, as fhhos4::Solver takes it (which checks it): the vectors of MatSetNearNullSpace(), or the constant one
	// (has_cnst, no vector: all the coefficients 1, the interpolation of 1 with one DoF per cell and face only)
	PetscErrorCode NearNullSpace(Mat P, IS cells, IS faces, Eigen::MatrixXd& cellNearNullSpace, Eigen::MatrixXd& faceNearNullSpace)
	{
		MPI_Comm comm = PetscObjectComm((PetscObject)P);
		MatNullSpace nearNullSpace;
		PetscBool hasConstant;
		PetscInt n;
		const Vec* vectors;
		PetscFunctionBegin;
		PetscCall(MatGetNearNullSpace(P, &nearNullSpace));
		PetscCheck(nearNullSpace, comm, PETSC_ERR_ARG_WRONGSTATE, "fhhos4: give the near null space of the operator, MatSetNearNullSpace(): the interpolation of 1 on the DoFs (a vector-valued unknown: of the constant of each component; elasticity: the rigid body modes). No default: a wrong one slows the AMG down or makes it diverge");
		PetscCall(MatNullSpaceGetVecs(nearNullSpace, &hasConstant, &n, &vectors));
		for (auto [is, matrix] : { std::pair{ cells, &cellNearNullSpace }, std::pair{ faces, &faceNearNullSpace } })
		{
			PetscInt size;
			PetscCall(ISGetLocalSize(is, &size));
			*matrix = Eigen::MatrixXd::Ones(size, n > 0 ? n : 1);
			for (PetscInt j = 0; j < n; j++)
			{
				Vec v;
				const PetscScalar* values;
				PetscCall(VecGetSubVector(vectors[j], is, &v));
				PetscCall(VecGetArrayRead(v, &values));
				matrix->col(j) = Eigen::Map<const Eigen::VectorXd>(values, size);
				PetscCall(VecRestoreArrayRead(v, &values));
				PetscCall(VecRestoreSubVector(vectors[j], is, &v));
			}
		}
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	// The options of the parameters in the messages of the inference of the discretization
	const DiscretizationInference::ParameterNames Options = { "-pc_fhhos4_dimension", "-pc_fhhos4_polynomial_space", "-pc_fhhos4_face_degree", "-pc_fhhos4_cell_degree",
	                                                          "-pc_fhhos4_vector_components", "p", "q" };

	// The parameters of the discretization whose options are not given (Dimension, FaceDegree, CellDegree -1,
	// PolynomialSpace::Inferred), inferred from the stored entries of the diagonal parts of A_TT and A_TF (the cells and
	// the own faces of each process; Solver/DiscretizationInference.h, ADR-0023) and written into the solver;
	// -pc_fhhos4_vector_components required unless the dimension and the degrees are given. Collective.
	PetscErrorCode InferDiscretization(fhhos4::Solver& s, Mat A_TT, Mat A_TF, MPI_Comm comm)
	{
		using namespace DiscretizationInference;
		Discretization set = { s.Dimension, s.PolynomialSpace, s.FaceDegree, s.CellDegree };
		bool degreesSet = s.Dimension != -1 && s.FaceDegree != -1 && s.CellDegree != -1;
		Block tt, tf;
		BlockAnalysis local;
		PetscMPIInt size;
		std::string error;
		bool allOk;
		PetscFunctionBegin;
		if (degreesSet && s.PolynomialSpace != fhhos4::PolynomialSpace::Inferred)
			PetscFunctionReturn(PETSC_SUCCESS);
		int components = s.VectorComponents >= 1 ? s.VectorComponents : 1;
		if (degreesSet) // the space alone: no analysis if P and Q have the same DoFs and solver (e.g. k = 0)
		{
			Inference ofP = DiscretizationInference::InferDiscretization(components * (int)PolynomialSpaceDimension(s.Dimension, s.CellDegree, false),
			                                                             components * (int)PolynomialSpaceDimension(s.Dimension - 1, s.FaceDegree, false), components, set, s.PCoarsening);
			if (ofP.Solvers == 1 && std::find(ofP.NotNeeded.begin(), ofP.NotNeeded.end(), Parameter::PolynomialSpace) != ofP.NotNeeded.end())
			{
				s.PolynomialSpace = fhhos4::PolynomialSpace::P;
				s.VectorComponents = components;
				PetscFunctionReturn(PETSC_SUCCESS);
			}
		}
		PetscCheck(degreesSet || s.VectorComponents >= 1, comm, PETSC_ERR_ARG_WRONGSTATE,
		           "fhhos4: give -pc_fhhos4_vector_components (1 for a scalar unknown): -pc_fhhos4_dimension, -pc_fhhos4_face_degree or -pc_fhhos4_cell_degree is not given, inferred from the blocks, whose DoFs per cell and face don't give the number of components");
		PetscCall(tt.Get(A_TT));
		PetscCall(tf.Get(A_TF));
		try
		{
			local = AnalyzeBlocks(tt.Diag.Eigen(), tf.Diag.Eigen());
		}
		catch (const fhhos4::Error& e)
		{
			error = e.what();
		}
		PetscCall(tt.Restore());
		PetscCall(tf.Restore());
		PetscCall(AllOk(comm, error.empty(), &allOk));
		if (!allOk && error.empty())
			error = "fhhos4: an input is wrong on another process (see its error)";
		Inference inference;
		if (allOk)
		{
			PetscCallMPI(MPI_Comm_size(comm, &size));
			std::vector<BlockAnalysis> analyses(size);
			PetscCallMPI(MPI_Allgather(&local, sizeof(BlockAnalysis), MPI_BYTE, analyses.data(), sizeof(BlockAnalysis), MPI_BYTE, comm));
			try // the same on all the processes
			{
				auto [cellBlockSize, faceBlockSize] = CombineAnalyses(analyses, Options);
				inference = DiscretizationInference::InferDiscretization(cellBlockSize, faceBlockSize, components, set, s.PCoarsening);
				if (inference.Solvers != 1)
					error = Message(inference, Options);
			}
			catch (const fhhos4::Error& e)
			{
				error = e.what();
			}
		}
		if (degreesSet && !error.empty()) // the stored entries don't give the space: P, as Solver::Setup() (on all the processes)
		{
			s.PolynomialSpace = fhhos4::PolynomialSpace::P;
			s.VectorComponents = components;
			PetscFunctionReturn(PETSC_SUCCESS);
		}
		PetscCheck(error.empty(), comm, PETSC_ERR_ARG_WRONG, "%s", error.c_str());
		s.Dimension = inference.Result.Dimension;
		s.PolynomialSpace = inference.Result.PolynomialSpace;
		s.FaceDegree = inference.Result.FaceDegree;
		s.CellDegree = inference.Result.CellDegree;
		if (s.VectorComponents == -1)
			s.VectorComponents = 1;
		if (s.Verbosity >= 1)
			PetscCall(PetscPrintf(comm, "%s\n", Description(inference, Options).c_str()));
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	PetscErrorCode PCSetUp_FHHOS4(PC pc)
	{
		PC_FHHOS4* ctx = (PC_FHHOS4*)pc->data;
		fhhos4::Solver& solver = ctx->Solver;
		MPI_Comm comm = PetscObjectComm((PetscObject)pc);
		Mat P = pc->pmat, A_TT, A_TF, A_FT = nullptr, A_FF;
		PetscMPIInt size;
		PetscBool nest, set, symmetric;
		PetscInt nRows, nCols, tBegin, tEnd, fBegin, fEnd;
		IS rows[2];
		int cbs = 0, fbs = 0;
		Block tt, tf, ft, ff;
		Eigen::MatrixXd cellNearNullSpace, faceNearNullSpace;
		std::string error;
		bool allOk;
		PetscFunctionBegin;
		PetscCallMPI(MPI_Comm_size(comm, &size));
#if FHHOS4_MPI
		solver.Comm = comm;
#else
		PetscCheck(size == 1, comm, PETSC_ERR_SUP, "fhhos4: the library is compiled without MPI (-DENABLE_MPI=OFF): one process only");
#endif
		PetscCall(PetscObjectTypeCompare((PetscObject)P, MATNEST, &nest));
		PetscCheck(nest, comm, PETSC_ERR_SUP, "fhhos4: the operator must be a MATNEST [[A_TT, A_TF], [A_FT, A_FF]] of AIJ blocks (the cells, then the faces), not %s", ((PetscObject)P)->type_name);
		PetscCall(MatNestGetSize(P, &nRows, &nCols));
		PetscCheck(nRows == 2 && nCols == 2, comm, PETSC_ERR_SUP, "fhhos4: the operator must be a MATNEST of 2 x 2 blocks [[A_TT, A_TF], [A_FT, A_FF]] (the cells, then the faces), not %" PetscInt_FMT " x %" PetscInt_FMT, nRows, nCols);
		symmetric = ctx->Symmetric;
		if (!ctx->SymmetricSet)
		{
			PetscCall(MatIsSymmetricKnown(P, &set, &symmetric));
			if (!set || !symmetric)
			{
				PetscCall(MatIsSPDKnown(P, &set, &symmetric));
				symmetric = (PetscBool)(set && symmetric);
			}
		}
		PetscCall(SubMatrix(P, 0, 0, "A_TT", &A_TT));
		PetscCall(SubMatrix(P, 0, 1, "A_TF", &A_TF));
		PetscCall(SubMatrix(P, 1, 1, "A_FF", &A_FF));
		if (!symmetric)
			PetscCall(SubMatrix(P, 1, 0, "A_FT", &A_FT));

		// The DoFs of each process: whole cells and faces
		PetscCall(InferDiscretization(solver, A_TT, A_TF, comm));
		try
		{
			cbs = solver.CellBlockSize();
			fbs = solver.FaceBlockSize();
		}
		catch (const fhhos4::Error& e)
		{
			SETERRQ(comm, PETSC_ERR_ARG_WRONGSTATE, "%s (options -pc_fhhos4_dimension, -pc_fhhos4_face_degree, -pc_fhhos4_cell_degree, -pc_fhhos4_vector_components)", e.what());
		}
		PetscCall(MatGetOwnershipRange(A_TT, &tBegin, &tEnd));
		PetscCall(MatGetOwnershipRange(A_FF, &fBegin, &fEnd));
		PetscCall(CheckOnAll(comm, tBegin % cbs == 0 && (tEnd - tBegin) % cbs == 0 && fBegin % fbs == 0 && (fEnd - fBegin) % fbs == 0,
		                     "each process must own the DoFs of whole cells and faces (" + std::to_string(cbs) + " DoFs per cell, " + std::to_string(fbs) + " per face): the first of its rows of A_TT and A_FF and their numbers multiples of these"));
		PetscCall(MatNestGetISs(P, rows, nullptr));
		PetscCall(PetscObjectReference((PetscObject)rows[0]));
		PetscCall(PetscObjectReference((PetscObject)rows[1]));
		PetscCall(ISDestroy(&ctx->Cells));
		PetscCall(ISDestroy(&ctx->Faces));
		ctx->Cells = rows[0];
		ctx->Faces = rows[1];
		PetscCall(NearNullSpace(P, ctx->Cells, ctx->Faces, cellNearNullSpace, faceNearNullSpace));

		// The blocks, read in place: on one process, the blocks themselves; on several, A_TT's diagonal part (cells are
		// not shared) and the others merged
		PetscCall(tt.Get(A_TT));
		PetscCall(tf.Get(A_TF));
		PetscCall(ff.Get(A_FF));
		if (A_FT)
			PetscCall(ft.Get(A_FT));
		if (size > 1)
		{
			try
			{
				if (tt.Offd.Eigen().nonZeros() > 0)
					error = "A_TT has entries in the columns of the cells of other processes: a cell belongs to one process";
				else if (symmetric)
					ctx->Merged = fhhos4::MergeOffDiag(tf.Diag.Eigen(), tf.Offd.Eigen(), tf.OffdColumns, ff.Diag.Eigen(), ff.Offd.Eigen(), ff.OffdColumns, fbs);
				else
					ctx->Merged = fhhos4::MergeOffDiag(tf.Diag.Eigen(), tf.Offd.Eigen(), tf.OffdColumns, ft.Diag.Eigen(), ft.Offd.Eigen(), ft.OffdColumns, ff.Diag.Eigen(), ff.Offd.Eigen(), ff.OffdColumns, fbs, cbs);
			}
			catch (const fhhos4::Error& e)
			{
				error = e.what();
			}
		}
		PetscCall(AllOk(comm, error.empty(), &allOk)); // the merge's errors are local
		if (allOk)
		{
			try
			{
				LibraryMatrix a_TT(tt.Diag.Eigen());
				if (size > 1 && symmetric)
					solver.SetupSym(a_TT, ctx->Merged.A_TF, ctx->Merged.A_FF, cellNearNullSpace, faceNearNullSpace, ctx->Merged.FaceGhostGlobalNumbers);
				else if (size > 1)
					solver.SetupNonSym(a_TT, ctx->Merged.A_TF, ctx->Merged.A_FT, ctx->Merged.A_FF, cellNearNullSpace, faceNearNullSpace, ctx->Merged.FaceGhostGlobalNumbers, ctx->Merged.CellGhostGlobalNumbers);
				else if (symmetric)
					solver.SetupSym(a_TT, LibraryMatrix(tf.Diag.Eigen()), LibraryMatrix(ff.Diag.Eigen()), cellNearNullSpace, faceNearNullSpace);
				else
				{
					ctx->A_FT = ft.Diag.Eigen();
					solver.SetupNonSym(a_TT, LibraryMatrix(tf.Diag.Eigen()), ctx->A_FT, LibraryMatrix(ff.Diag.Eigen()), cellNearNullSpace, faceNearNullSpace);
				}
				ctx->SetUpSymmetric = symmetric;
			}
			catch (const fhhos4::Error& e) // collective: raised on all the processes
			{
				error = e.what();
			}
		}
		for (Block* block : { &tt, &tf, &ft, &ff })
			PetscCall(block->Restore());
		PetscCheck(allOk, comm, PETSC_ERR_ARG_WRONG, "fhhos4: %s", error.empty() ? "an input is wrong on another process (see its error)" : error.c_str());
		PetscCheck(error.empty(), comm, PETSC_ERR_LIB, "fhhos4: %s", error.c_str());
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	// z = M^{-1} r: the cells and faces of r and z (VecGetSubVector() with the index sets of the operator's rows)
	PetscErrorCode PCApply_FHHOS4(PC pc, Vec r, Vec z)
	{
		PC_FHHOS4* ctx = (PC_FHHOS4*)pc->data;
		Vec r_T, r_F, z_T, z_F;
		const PetscScalar *rT, *rF;
		PetscScalar *zT, *zF;
		PetscInt nT, nF;
		std::string error;
		PetscFunctionBegin;
		PetscCall(VecGetSubVector(r, ctx->Cells, &r_T));
		PetscCall(VecGetSubVector(r, ctx->Faces, &r_F));
		PetscCall(VecGetSubVector(z, ctx->Cells, &z_T));
		PetscCall(VecGetSubVector(z, ctx->Faces, &z_F));
		PetscCall(VecGetLocalSize(r_T, &nT));
		PetscCall(VecGetLocalSize(r_F, &nF));
		PetscCall(VecGetArrayRead(r_T, &rT));
		PetscCall(VecGetArrayRead(r_F, &rF));
		PetscCall(VecGetArrayWrite(z_T, &zT));
		PetscCall(VecGetArrayWrite(z_F, &zF));
		{
			Eigen::Map<const Eigen::VectorXd> b_T(rT, nT), b_F(rF, nF);
			Eigen::Map<Eigen::VectorXd> x_T(zT, nT), x_F(zF, nF);
			x_F.setZero(); // the initial guess of the cycle
			try
			{
				ctx->Solver.Solve(b_T, b_F, x_T, x_F);
			}
			catch (const fhhos4::Error& e)
			{
				error = e.what();
			}
		}
		PetscCall(VecRestoreArrayWrite(z_F, &zF));
		PetscCall(VecRestoreArrayWrite(z_T, &zT));
		PetscCall(VecRestoreArrayRead(r_F, &rF));
		PetscCall(VecRestoreArrayRead(r_T, &rT));
		PetscCall(VecRestoreSubVector(z, ctx->Faces, &z_F));
		PetscCall(VecRestoreSubVector(z, ctx->Cells, &z_T));
		PetscCall(VecRestoreSubVector(r, ctx->Faces, &r_F));
		PetscCall(VecRestoreSubVector(r, ctx->Cells, &r_T));
		PetscCheck(error.empty(), PetscObjectComm((PetscObject)pc), PETSC_ERR_LIB, "fhhos4: %s", error.c_str());
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	const char* const Cycles[] = { "v", "w", "k" };      // fhhos4::Cycle
	const char* const Krylovs[] = { "none", "fcg", "gcr" }; // fhhos4::Krylov None, FCG, GCR
	const char* const Spaces[] = { "p", "q", "inferred" }; // fhhos4::PolynomialSpace

	PetscErrorCode PCSetFromOptions_FHHOS4(PC pc, PetscOptionItems PetscOptionsObject)
	{
		PC_FHHOS4* ctx = (PC_FHHOS4*)pc->data;
		fhhos4::Solver& s = ctx->Solver;
		PetscInt value, krylov = s.Krylov == fhhos4::Krylov::FCG ? 1 : s.Krylov == fhhos4::Krylov::GCR ? 2 : 0;
		PetscReal tolerance;
		PetscBool set;
		PetscFunctionBegin;
		PetscOptionsHeadBegin(PetscOptionsObject, "fhhos4 options (U-AMG for hybrid discretizations)");
		PetscCall(PetscOptionsInt("-pc_fhhos4_dimension", "Dimension of the mesh: 2 or 3 (default: inferred from the blocks)", "PCFHHOS4", s.Dimension, &value, &set));
		if (set)
			s.Dimension = (int)value;
		PetscCall(PetscOptionsInt("-pc_fhhos4_face_degree", "Degree of the polynomials on the faces (default: inferred from the blocks)", "PCFHHOS4", s.FaceDegree, &value, &set));
		if (set)
			s.FaceDegree = (int)value;
		PetscCall(PetscOptionsInt("-pc_fhhos4_cell_degree", "Degree of the polynomials in the cells: k - 1, k or k + 1 (default: inferred from the blocks)", "PCFHHOS4", s.CellDegree, &value, &set));
		if (set)
			s.CellDegree = (int)value;
		PetscCall(PetscOptionsEList("-pc_fhhos4_polynomial_space", "Space of the cell and face bases: p (degree <= k), q (degree <= k in each variable), inferred from the blocks (the default)", "PCFHHOS4", Spaces, 3, Spaces[(int)s.PolynomialSpace], &value, &set));
		if (set)
			s.PolynomialSpace = (fhhos4::PolynomialSpace)value;
		PetscCall(PetscOptionsInt("-pc_fhhos4_vector_components", "Components of the unknown, interleaved in each block (required if the dimension or a degree is inferred; default: 1 otherwise)", "PCFHHOS4", s.VectorComponents, &value, &set));
		if (set)
			s.VectorComponents = (int)value;
		PetscCall(PetscOptionsBool("-pc_fhhos4_symmetric", "Symmetric matrix (default: MatIsSymmetricKnown() or MatIsSPDKnown() of the operator)", "PCFHHOS4", ctx->Symmetric, &ctx->Symmetric, &ctx->SymmetricSet));
		PetscCall(PetscOptionsEList("-pc_fhhos4_cycle", "Cycle: k (non-linear: a flexible Krylov method, fcg, fgmres or gcr), v, w", "PCFHHOS4", Cycles, 3, Cycles[(int)s.Cycle], &value, &set));
		if (set)
			s.Cycle = (fhhos4::Cycle)value;
		PetscCall(PetscOptionsEList("-pc_fhhos4_krylov", "fhhos4's Krylov method in each application: none (one cycle), fcg, gcr (with -ksp_type preonly)", "PCFHHOS4", Krylovs, 3, Krylovs[krylov], &value, &set));
		if (set)
			s.Krylov = value == 1 ? fhhos4::Krylov::FCG : value == 2 ? fhhos4::Krylov::GCR : fhhos4::Krylov::None;
		PetscCall(PetscOptionsInt("-pc_fhhos4_max_it", "Iterations of fhhos4's Krylov method (cycles with -pc_fhhos4_krylov none)", "PCFHHOS4", s.MaxIterations, &value, &set));
		if (set)
			s.MaxIterations = (int)value;
		PetscCall(PetscOptionsReal("-pc_fhhos4_rtol", "Relative tolerance of fhhos4's Krylov method", "PCFHHOS4", s.Tolerance, &tolerance, &set));
		if (set)
			s.Tolerance = tolerance;
		PetscCall(PetscOptionsInt("-pc_fhhos4_verbosity", "fhhos4's printouts: 0, 1 (one line per setup and solve), 2", "PCFHHOS4", s.Verbosity, &value, &set));
		if (set)
			s.Verbosity = (int)value;
		PetscOptionsHeadEnd();
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	PetscErrorCode PCView_FHHOS4(PC pc, PetscViewer viewer)
	{
		PC_FHHOS4* ctx = (PC_FHHOS4*)pc->data;
		const fhhos4::Solver& s = ctx->Solver;
		PetscBool ascii;
		PetscFunctionBegin;
		PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &ascii));
		if (ascii)
		{
			PetscCall(PetscViewerASCIIPrintf(viewer, "  fhhos4 U-AMG: Dimension %d, FaceDegree %d, CellDegree %d, PolynomialSpace %s, VectorComponents %d\n", s.Dimension, s.FaceDegree, s.CellDegree,
			                                 s.PolynomialSpace == fhhos4::PolynomialSpace::Q ? "Q" : s.PolynomialSpace == fhhos4::PolynomialSpace::P ? "P" : "inferred", s.VectorComponents));
			PetscCall(PetscViewerASCIIPrintf(viewer, "  %s-cycle%s\n", s.Cycle == fhhos4::Cycle::K ? "K" : s.Cycle == fhhos4::Cycle::V ? "V" : "W",
			                                 s.Cycle == fhhos4::Cycle::K ? ": non-linear, the Krylov method must be flexible (fcg, fgmres, gcr)" : ""));
			if (s.Krylov == fhhos4::Krylov::None)
				PetscCall(PetscViewerASCIIPrintf(viewer, "  %d cycle(s) per application\n", s.MaxIterations));
			else
				PetscCall(PetscViewerASCIIPrintf(viewer, "  fhhos4's %s in each application: %d iterations at most, relative tolerance %g\n", s.Krylov == fhhos4::Krylov::FCG ? "FCG" : "GCR", s.MaxIterations, s.Tolerance));
			if (s.IsSetUp())
				PetscCall(PetscViewerASCIIPrintf(viewer, "  %s setup, %d levels\n", ctx->SetUpSymmetric ? "symmetric" : "non-symmetric", s.NumberOfLevels()));
		}
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	PetscErrorCode PCReset_FHHOS4(PC pc)
	{
		PC_FHHOS4* ctx = (PC_FHHOS4*)pc->data;
		PetscFunctionBegin;
		PetscCall(ISDestroy(&ctx->Cells));
		PetscCall(ISDestroy(&ctx->Faces));
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	PetscErrorCode PCDestroy_FHHOS4(PC pc)
	{
		PetscFunctionBegin;
		PetscCall(PCReset_FHHOS4(pc));
		delete (PC_FHHOS4*)pc->data;
		pc->data = nullptr;
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	PetscErrorCode PCCreate_FHHOS4(PC pc)
	{
		PC_FHHOS4* ctx = new (std::nothrow) PC_FHHOS4();
		PetscFunctionBegin;
		PetscCheck(ctx, PetscObjectComm((PetscObject)pc), PETSC_ERR_MEM, "fhhos4: out of memory");
		ctx->Solver.Krylov = fhhos4::Krylov::None; // one cycle per application
		ctx->Solver.MaxIterations = 1;
		pc->data = ctx;
		pc->ops->setup = PCSetUp_FHHOS4;
		pc->ops->apply = PCApply_FHHOS4;
		pc->ops->setfromoptions = PCSetFromOptions_FHHOS4;
		pc->ops->view = PCView_FHHOS4;
		pc->ops->reset = PCReset_FHHOS4;
		pc->ops->destroy = PCDestroy_FHHOS4;
		PetscFunctionReturn(PETSC_SUCCESS);
	}
}

PETSC_EXTERN PetscErrorCode PetscDLLibraryRegister_fhhos4_petsc(void)
{
	PetscFunctionBegin;
	PetscCall(PCRegister(PCFHHOS4, PCCreate_FHHOS4));
	PetscFunctionReturn(PETSC_SUCCESS);
}
