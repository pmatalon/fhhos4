// fhhos4's preconditioner for PETSc (library/petsc, -pc_type fhhos4), loaded by -dll_append as by any PETSc program
// (ctest PETScPC_np1...4): the Cartesian hybrid system partitioned among the processes (PartitionedSystem), each
// process giving the rows of its cells and faces to a MATNEST of AIJ blocks (MATSEQAIJ on one process, MATMPIAIJ on
// several), solved by PETSc's Krylov methods preconditioned by fhhos4, against the library's solve of the whole system
// on one process; or as one AIJ matrix with the index sets of its cells and faces (PCFHHOS4SetIS()). The results
// printed by the first process; a test fails if it fails on one process.
#include <algorithm>
#include <numeric>
#include <random>
#include <sstream>
#include <gtest/gtest.h>
#include <petscksp.h>
#include <fhhos4>
#include "fhhos4_petsc.h"
#include "support/CartesianHybridSystem.h"
#include "support/Kronecker.h"
#include "support/PartitionedSystem.h"

using namespace fhhos4_tests;

// The near null space of c components at degree 0 (c DoFs per cell and face): the constant of each component, one
// column each (1 on its DoF of each block)
static Eigen::MatrixXd Constants(Eigen::Index nBlocks, int c)
{
	return Eigen::MatrixXd::Identity(c, c).replicate(nBlocks, 1);
}

namespace
{
	const Communicator& World()
	{
		return Communicator::Current();
	}

	// The vectors of all the processes, in rank order
	Vector Gathered(const Vector& local)
	{
		vector<double> values(local.data(), local.data() + local.rows());
		vector<double> all = World().AllGatherV(values);
		return Eigen::Map<Vector>(all.data(), all.size());
	}

	// An AIJ matrix from the rows of the process with their global columns
	Mat AIJ(const SparseMatrix& rows, PetscInt nOwnCols)
	{
		Mat A;
		vector<PetscInt> i(rows.outerIndexPtr(), rows.outerIndexPtr() + rows.rows() + 1), j(rows.innerIndexPtr(), rows.innerIndexPtr() + rows.nonZeros());
		PetscCallAbort(PETSC_COMM_WORLD, MatCreate(PETSC_COMM_WORLD, &A));
		PetscCallAbort(PETSC_COMM_WORLD, MatSetSizes(A, rows.rows(), nOwnCols, PETSC_DETERMINE, PETSC_DETERMINE));
		PetscCallAbort(PETSC_COMM_WORLD, MatSetType(A, MATAIJ));
		PetscCallAbort(PETSC_COMM_WORLD, MatSeqAIJSetPreallocationCSR(A, i.data(), j.data(), rows.valuePtr()));
		PetscCallAbort(PETSC_COMM_WORLD, MatMPIAIJSetPreallocationCSR(A, i.data(), j.data(), rows.valuePtr()));
		return A;
	}

	// The MATNEST of the rows of the process of the blocks of a global system (cellDoFs DoFs per cell, faceDoFs per face):
	// [[A_TT, A_TF], [A_FT, A_FF]], its AIJ blocks in blocks
	Mat Nest(const SparseMatrix (&globals)[4], const IndexMap& cells, const IndexMap& faces, int cellDoFs, int faceDoFs, Mat (&blocks)[4])
	{
		const IndexMap* maps[2] = { &cells, &faces };
		int dofs[2] = { cellDoFs, faceDoFs };
		for (int k = 0; k < 4; k++)
		{
			int row = k / 2, col = k % 2;
			blocks[k] = AIJ(globals[k].middleRows(dofs[row] * maps[row]->Begin, dofs[row] * maps[row]->NOwned), dofs[col] * maps[col]->NOwned);
		}
		Mat P;
		PetscCallAbort(PETSC_COMM_WORLD, MatCreateNest(PETSC_COMM_WORLD, 2, nullptr, 2, nullptr, blocks, &P));
		return P;
	}

	// The message of the last error raised by PETSc (StoreError(): returned, not printed)
	std::string LastError;

	PetscErrorCode StoreError(MPI_Comm, int, const char*, const char*, PetscErrorCode n, PetscErrorType p, const char* message, void*)
	{
		if (p == PETSC_ERROR_INITIAL)
			LastError = message;
		return n;
	}

	// The parts of a vector of the operator on the cells and faces of the process (rows: theirs in the operator)
	void SetParts(IS const (&rows)[2], Vec v, const Vector& v_T, const Vector& v_F)
	{
		for (auto [is, values] : { pair{ rows[0], &v_T }, pair{ rows[1], &v_F } })
		{
			Vec part;
			PetscScalar* a;
			PetscCallAbort(PETSC_COMM_WORLD, VecGetSubVector(v, is, &part));
			PetscCallAbort(PETSC_COMM_WORLD, VecGetArrayWrite(part, &a));
			Eigen::Map<Vector>(a, values->rows()) = *values;
			PetscCallAbort(PETSC_COMM_WORLD, VecRestoreArrayWrite(part, &a));
			PetscCallAbort(PETSC_COMM_WORLD, VecRestoreSubVector(v, is, &part));
		}
	}

	Vector Part(IS const (&rows)[2], Vec v, int i)
	{
		Vec part;
		const PetscScalar* a;
		PetscInt n;
		PetscCallAbort(PETSC_COMM_WORLD, VecGetSubVector(v, rows[i], &part));
		PetscCallAbort(PETSC_COMM_WORLD, VecGetLocalSize(part, &n));
		PetscCallAbort(PETSC_COMM_WORLD, VecGetArrayRead(part, &a));
		Vector values = Eigen::Map<const Vector>(a, n);
		PetscCallAbort(PETSC_COMM_WORLD, VecRestoreArrayRead(part, &a));
		PetscCallAbort(PETSC_COMM_WORLD, VecRestoreSubVector(v, rows[i], &part));
		return values;
	}

	enum class NearNullSpace { InterpolationOfOne, Constant, None };

	// The system of n x n cells (non-symmetric if asymmetry != 0, on c components: the blocks ⊗ I_c) as a MATNEST of the
	// rows of each process, with the right-hand side, the interpolation of 1 as near null space (PETSc's orthonormal
	// vectors: one per component, normalized), and the solution and iterations of the library on one process. The
	// vectors: those of MatCreateVecs() (plain, with the layout of the nest), or nestVectors, VECNEST ones (a sub-vector
	// for the cells, one for the faces). MakeMonolithic(): the operator as one AIJ matrix instead.
	struct PETScSystem
	{
		int c;
		Mat Blocks[4] = {}, P = nullptr;
		IS Rows[2] = {};            // of the cells and of the faces of the process in P
		IS GivenFaces = nullptr;    // P one AIJ matrix: the faces given to the PC (as PCFHHOS4SetIS()), Rows[1] or fewer
		Eigen::Index CellOffset, FaceOffset; // of the DoFs of the process in the global vectors
		Vec b = nullptr;
		Vector ExpectedX_T, ExpectedX_F;
		int ExpectedIterations = 0;       // by the setup of the system's symmetry
		int ExpectedNonSymIterations = 0; // by the non-symmetric setup

		PETScSystem(double asymmetry, int components, bool symmetricFlag, NearNullSpace nearNullSpace = NearNullSpace::InterpolationOfOne, bool nestVectors = false)
			: c(components)
		{
			int n = 60;
			CartesianHybridSystem global(n, asymmetry);
			PartitionedSystem system(global, n, PartitionShape::Blocks, World());
			const IndexMap& cells = *system.Cells;
			const IndexMap& faces = *system.Faces;
			Vector b_T = OnEachComponent(Vector::LinSpaced(global.A_TT.rows(), 1, 2), c), b_F = Vector::Zero(c * global.S.rows());
			Eigen::MatrixXd cellOnes = Constants(global.A_TT.rows(), c), faceOnes = Constants(global.S.rows(), c);
			fhhos4::SparseMatrix A_TT = OnEachComponent(system.Global.A_TT, c), A_TF = OnEachComponent(system.Global.A_TF, c);
			fhhos4::SparseMatrix A_FT = OnEachComponent(system.Global.A_FT, c), A_FF = OnEachComponent(system.Global.A_FF, c);
			auto reference = [&](bool symmetric, Vector& x_T, Vector& x_F)
			{
				fhhos4::Solver solver;
				solver.Dimension = 2;
				solver.FaceDegree = 0;
				solver.CellDegree = 0;
				solver.Tolerance = 1e-10;
				solver.VectorComponents = c;
				if (symmetric)
					solver.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes);
				else
					solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellOnes, faceOnes);
				x_T.resize(b_T.rows());
				x_F = Vector::Zero(b_F.rows());
				return solver.Solve(b_T, b_F, x_T, x_F).Iterations;
			};
			ExpectedIterations = reference(asymmetry == 0, ExpectedX_T, ExpectedX_F);
			Vector x_T, x_F;
			ExpectedNonSymIterations = asymmetry == 0 ? reference(false, x_T, x_F) : ExpectedIterations;
			SparseMatrix globals[4] = { OnEachComponent(system.Global.A_TT, c), OnEachComponent(system.Global.A_TF, c), OnEachComponent(system.Global.A_FT, c), OnEachComponent(system.Global.A_FF, c) };
			P = Nest(globals, cells, faces, c, c, Blocks);
			PetscCallAbort(PETSC_COMM_WORLD, MatNestGetISs(P, Rows, nullptr));
			for (IS rows : Rows)
				PetscCallAbort(PETSC_COMM_WORLD, PetscObjectReference((PetscObject)rows));
			CellOffset = c * cells.Begin;
			FaceOffset = c * faces.Begin;
			if (symmetricFlag)
				PetscCallAbort(PETSC_COMM_WORLD, MatSetOption(P, MAT_SYMMETRIC, PETSC_TRUE));
			if (nestVectors)
			{
				Vec parts[2];
				PetscCallAbort(PETSC_COMM_WORLD, MatCreateVecs(Blocks[0], nullptr, &parts[0]));
				PetscCallAbort(PETSC_COMM_WORLD, MatCreateVecs(Blocks[3], nullptr, &parts[1]));
				PetscCallAbort(PETSC_COMM_WORLD, VecCreateNest(PETSC_COMM_WORLD, 2, nullptr, parts, &b));
				for (Vec& v : parts)
					PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&v));
			}
			else
				PetscCallAbort(PETSC_COMM_WORLD, MatCreateVecs(P, nullptr, &b));
			SetParts(Rows, b, b_T.segment(CellOffset, c * cells.NOwned), b_F.segment(FaceOffset, c * faces.NOwned));
			if (nearNullSpace != NearNullSpace::None)
			{
				vector<Vec> ones(nearNullSpace == NearNullSpace::Constant ? 0 : c);
				for (int k = 0; k < (int)ones.size(); k++)
				{
					// The component k of 1: the first function of each block on that component
					auto onComponent = [&](BigNumber blocks) { Vector v = Vector::Zero(c * blocks); v(Eigen::seqN(k, blocks, c)).setOnes(); return v; };
					PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(b, &ones[k]));
					SetParts(Rows, ones[k], onComponent(cells.NOwned), onComponent(faces.NOwned));
					PetscCallAbort(PETSC_COMM_WORLD, VecNormalize(ones[k], nullptr));
				}
				MatNullSpace space;
				PetscCallAbort(PETSC_COMM_WORLD, MatNullSpaceCreate(PETSC_COMM_WORLD, nearNullSpace == NearNullSpace::Constant ? PETSC_TRUE : PETSC_FALSE, ones.size(), ones.data(), &space));
				PetscCallAbort(PETSC_COMM_WORLD, MatSetNearNullSpace(P, space));
				PetscCallAbort(PETSC_COMM_WORLD, MatNullSpaceDestroy(&space));
				for (Vec& v : ones)
					PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&v));
			}
		}

		~PETScSystem()
		{
			for (IS& rows : Rows)
				PetscCallAbort(PETSC_COMM_WORLD, ISDestroy(&rows));
			PetscCallAbort(PETSC_COMM_WORLD, ISDestroy(&GivenFaces));
			PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&b));
			PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&P));
			for (Mat& A : Blocks)
				PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&A));
		}

		// The operator as one AIJ matrix (MatConvert() of the MATNEST), its rows shuffled on each process (the DoFs of the
		// cells and faces mixed, Rows following them), given to the PC with PCFHHOS4SetIS(); the right-hand side and the
		// near null space permuted alike (VecPermute()). With leftOutEvery > 0, every leftOutEvery-th face of the process
		// eliminated as Dirichlet DoFs, MatZeroRowsColumns() with the expected solution (the system's solution unchanged),
		// and left out of GivenFaces.
		void MakeMonolithic(int leftOutEvery)
		{
			Mat A, permuted;
			IS permutation;
			PetscInt begin, end;
			PetscBool set, symmetric;
			MatNullSpace nearNullSpace;
			PetscCallAbort(PETSC_COMM_WORLD, MatConvert(P, MATAIJ, MAT_INITIAL_MATRIX, &A));
			PetscCallAbort(PETSC_COMM_WORLD, MatGetOwnershipRange(A, &begin, &end));
			vector<PetscInt> order(end - begin), position(end - begin);
			std::iota(order.begin(), order.end(), begin);
			std::shuffle(order.begin(), order.end(), std::mt19937(World().Rank()));
			for (PetscInt k = 0; k < end - begin; k++)
				position[order[k] - begin] = begin + k;
			PetscCallAbort(PETSC_COMM_WORLD, ISCreateGeneral(PETSC_COMM_WORLD, end - begin, order.data(), PETSC_COPY_VALUES, &permutation));
			PetscCallAbort(PETSC_COMM_WORLD, MatPermute(A, permutation, permutation, &permuted)); // its row k: the row order[k] of A
			for (IS& rows : Rows)
			{
				const PetscInt* indices;
				PetscInt n;
				PetscCallAbort(PETSC_COMM_WORLD, ISGetLocalSize(rows, &n));
				PetscCallAbort(PETSC_COMM_WORLD, ISGetIndices(rows, &indices));
				vector<PetscInt> moved(n);
				for (PetscInt k = 0; k < n; k++)
					moved[k] = position[indices[k] - begin];
				PetscCallAbort(PETSC_COMM_WORLD, ISRestoreIndices(rows, &indices));
				PetscCallAbort(PETSC_COMM_WORLD, ISDestroy(&rows));
				PetscCallAbort(PETSC_COMM_WORLD, ISCreateGeneral(PETSC_COMM_WORLD, n, moved.data(), PETSC_COPY_VALUES, &rows));
			}
			PetscCallAbort(PETSC_COMM_WORLD, VecPermute(b, permutation, PETSC_FALSE));
			PetscCallAbort(PETSC_COMM_WORLD, MatGetNearNullSpace(P, &nearNullSpace));
			if (nearNullSpace)
			{
				PetscBool hasConstant;
				PetscInt n;
				const Vec* vectors;
				MatNullSpace permutedSpace;
				PetscCallAbort(PETSC_COMM_WORLD, MatNullSpaceGetVecs(nearNullSpace, &hasConstant, &n, &vectors));
				vector<Vec> moved(n);
				for (PetscInt k = 0; k < n; k++)
				{
					PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(vectors[k], &moved[k]));
					PetscCallAbort(PETSC_COMM_WORLD, VecCopy(vectors[k], moved[k]));
					PetscCallAbort(PETSC_COMM_WORLD, VecPermute(moved[k], permutation, PETSC_FALSE));
				}
				PetscCallAbort(PETSC_COMM_WORLD, MatNullSpaceCreate(PETSC_COMM_WORLD, hasConstant, n, moved.data(), &permutedSpace));
				PetscCallAbort(PETSC_COMM_WORLD, MatSetNearNullSpace(permuted, permutedSpace));
				PetscCallAbort(PETSC_COMM_WORLD, MatNullSpaceDestroy(&permutedSpace));
				for (Vec& v : moved)
					PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&v));
			}
			if (leftOutEvery > 0)
			{
				Vec x;
				const PetscInt* faceRows;
				PetscInt nF, nT;
				vector<PetscInt> eliminated, kept;
				PetscCallAbort(PETSC_COMM_WORLD, ISGetLocalSize(Rows[0], &nT));
				PetscCallAbort(PETSC_COMM_WORLD, ISGetLocalSize(Rows[1], &nF));
				PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(b, &x));
				SetParts(Rows, x, ExpectedX_T.segment(CellOffset, nT), ExpectedX_F.segment(FaceOffset, nF));
				PetscCallAbort(PETSC_COMM_WORLD, ISGetIndices(Rows[1], &faceRows));
				for (PetscInt f = 0; f < nF / c; f++)
				{
					vector<PetscInt>& to = f % leftOutEvery == 0 ? eliminated : kept;
					to.insert(to.end(), faceRows + c * f, faceRows + c * (f + 1));
				}
				PetscCallAbort(PETSC_COMM_WORLD, ISRestoreIndices(Rows[1], &faceRows));
				PetscCallAbort(PETSC_COMM_WORLD, MatZeroRowsColumns(permuted, eliminated.size(), eliminated.data(), 1.0, x, b));
				PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&x));
				PetscCallAbort(PETSC_COMM_WORLD, ISCreateGeneral(PETSC_COMM_WORLD, kept.size(), kept.data(), PETSC_COPY_VALUES, &GivenFaces));
			}
			else
			{
				GivenFaces = Rows[1];
				PetscCallAbort(PETSC_COMM_WORLD, PetscObjectReference((PetscObject)GivenFaces));
			}
			PetscCallAbort(PETSC_COMM_WORLD, MatIsSymmetricKnown(P, &set, &symmetric));
			if (set && symmetric)
				PetscCallAbort(PETSC_COMM_WORLD, MatSetOption(permuted, MAT_SYMMETRIC, PETSC_TRUE));
			PetscCallAbort(PETSC_COMM_WORLD, ISDestroy(&permutation));
			PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&A));
			PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&P));
			P = permuted;
		}

		struct Result
		{
			PetscErrorCode ErrorCode; // of KSPSolve(), with the setup
			PetscInt Iterations;
			KSPConvergedReason Reason;
			double Error; // relative, of the solution against the library's on one process
		};

		// The solve by the Krylov method kspType preconditioned by fhhos4 (-pc_type fhhos4), with the options for this
		// solve (prefix), the relative tolerance 1e-10 on the unpreconditioned residual (with a Krylov method)
		Result Solve(const char* kspType, const string& prefix, const string& options)
		{
			KSP ksp;
			PC pc;
			Vec x;
			Result result{};
			PetscCallAbort(PETSC_COMM_WORLD, KSPCreate(PETSC_COMM_WORLD, &ksp));
			PetscCallAbort(PETSC_COMM_WORLD, KSPSetOptionsPrefix(ksp, prefix.c_str()));
			PetscCallAbort(PETSC_COMM_WORLD, KSPSetOperators(ksp, P, P));
			PetscCallAbort(PETSC_COMM_WORLD, KSPSetType(ksp, kspType));
			PetscCallAbort(PETSC_COMM_WORLD, KSPGetPC(ksp, &pc));
			PetscCallAbort(PETSC_COMM_WORLD, PCSetType(pc, PCFHHOS4)); // registered by -dll_append
			if (GivenFaces) // PCFHHOS4SetIS() composes them so; the library not linked here (-dll_append)
			{
				PetscCallAbort(PETSC_COMM_WORLD, PetscObjectCompose((PetscObject)pc, "fhhos4_cells", (PetscObject)Rows[0]));
				PetscCallAbort(PETSC_COMM_WORLD, PetscObjectCompose((PetscObject)pc, "fhhos4_faces", (PetscObject)GivenFaces));
			}
			if (string(kspType) != KSPPREONLY)
			{
				PetscCallAbort(PETSC_COMM_WORLD, KSPSetTolerances(ksp, 1e-10, PETSC_CURRENT, PETSC_CURRENT, 200));
				PetscCallAbort(PETSC_COMM_WORLD, KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
			}
			istringstream in(options);
			string token, prefixed;
			while (in >> token)
				prefixed += (token[0] == '-' ? "-" + prefix + token.substr(1) : token) + " ";
			PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsInsertString(nullptr, prefixed.c_str()));
			PetscCallAbort(PETSC_COMM_WORLD, KSPSetFromOptions(ksp));
			PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(b, &x));
			// The setup in KSPSolve(), its work vectors like b: KSPSetUp() before would create them from the MATNEST (plain
			// vectors), into which PETSc can't copy VECNEST ones
			PetscCallAbort(PETSC_COMM_WORLD, PetscPushErrorHandler(StoreError, nullptr));
			result.ErrorCode = KSPSolve(ksp, b, x);
			PetscCallAbort(PETSC_COMM_WORLD, PetscPopErrorHandler());
			if (result.ErrorCode == PETSC_SUCCESS)
			{
				PetscCallAbort(PETSC_COMM_WORLD, KSPGetIterationNumber(ksp, &result.Iterations));
				PetscCallAbort(PETSC_COMM_WORLD, KSPGetConvergedReason(ksp, &result.Reason));
				Vector allX(ExpectedX_T.rows() + ExpectedX_F.rows()), expected(allX.rows());
				allX << Gathered(Part(Rows, x, 0)), Gathered(Part(Rows, x, 1));
				expected << ExpectedX_T, ExpectedX_F;
				result.Error = (allX - expected).norm() / expected.norm();
				if (World().IsRoot())
					cout << "\t" << World().Size() << " processes, " << kspType << " + fhhos4 (" << options << "): " << result.Iterations << " iterations, relative error " << result.Error << endl;
			}
			PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&x));
			PetscCallAbort(PETSC_COMM_WORLD, KSPDestroy(&ksp));
			return result;
		}
	};

	const string Degrees = " -pc_fhhos4_dimension 2 -pc_fhhos4_face_degree 0 -pc_fhhos4_cell_degree 0";

	// Converged to the solution, in the iterations of the library alone on one process (its Krylov method): within 1 on
	// one process (PETSc's FGMRES against the library's GCR, its residual norms against the condensed system's); within
	// 4 on several, as in DistributedLibrary (each process aggregates its own cells; the hybrid smoothers depend on the
	// partition, and on the threads of the processes, which share the cores)
	void ExpectSolves(const PETScSystem::Result& result, int expectedIterations)
	{
		if (World().IsRoot())
			cout << "\t\tthe library alone on one process: " << expectedIterations << " iterations" << endl;
		EXPECT_EQ(result.ErrorCode, PETSC_SUCCESS) << LastError;
		EXPECT_GT(result.Reason, 0);
		EXPECT_LE(result.Error, 1e-8);
		EXPECT_LE(abs(result.Iterations - expectedIterations), World().Size() == 1 ? 1 : 4);
	}
}

// A symmetric system (MAT_SYMMETRIC: the symmetric setup), PETSc's FCG with one K-cycle per application: the iterations
// of the library alone (its FCG with the K-cycle)
TEST(PETScPC, Symmetric)
{
	PETScSystem system(0, 1, true);
	ExpectSolves(system.Solve(KSPFCG, "sym_", Degrees), system.ExpectedIterations);
}

// No symmetry flag: the non-symmetric setup (with FGMRES), also of a symmetric system (the iterations of the library's
// non-symmetric setup, its l1 hybrid Gauss-Seidel and GCR, not those of the symmetric one); -pc_fhhos4_symmetric: the
// symmetric one
TEST(PETScPC, NonSymmetric)
{
	PETScSystem nonsymmetric(0.3, 1, false);
	ExpectSolves(nonsymmetric.Solve(KSPFGMRES, "nonsym_", Degrees), nonsymmetric.ExpectedIterations);
	PETScSystem symmetric(0, 1, false);
	ExpectSolves(symmetric.Solve(KSPFGMRES, "nonsymsym_", Degrees), symmetric.ExpectedNonSymIterations);
	ExpectSolves(symmetric.Solve(KSPFCG, "forcedsym_", Degrees + " -pc_fhhos4_symmetric"), symmetric.ExpectedIterations);
}

// A vector-valued unknown (-pc_fhhos4_vector_components 2: blocks of 2 DoFs, interleaved), the near null space one vector
// per component
TEST(PETScPC, VectorComponents)
{
	PETScSystem system(0, 2, true);
	ExpectSolves(system.Solve(KSPFCG, "comp_", Degrees + " -pc_fhhos4_vector_components 2"), system.ExpectedIterations);
}

// Nested vectors (VECNEST, a sub-vector for the cells, one for the faces) for the right-hand side, the solution and the
// near null space (2 components: two vectors), as plain ones (MatCreateVecs() of a MATNEST): the same iterations
TEST(PETScPC, NestVectors)
{
	PETScSystem plain(0, 2, true), nest(0, 2, true, NearNullSpace::InterpolationOfOne, true);
	PETScSystem::Result plainResult = plain.Solve(KSPFCG, "plain_", Degrees + " -pc_fhhos4_vector_components 2");
	PETScSystem::Result nestResult = nest.Solve(KSPFCG, "nestvectors_", Degrees + " -pc_fhhos4_vector_components 2");
	ExpectSolves(nestResult, nest.ExpectedIterations);
	EXPECT_EQ(nestResult.Iterations, plainResult.Iterations);
}

// -ksp_type preonly with fhhos4's FCG inside the preconditioner: the library's solve; PETSc's constant near null space
// (has_cnst), the interpolation of 1 at degree 0 on one component
TEST(PETScPC, PreonlyAndConstantNearNullSpace)
{
	PETScSystem system(0, 1, true, NearNullSpace::Constant);
	PETScSystem::Result result = system.Solve(KSPPREONLY, "preonly_", Degrees + " -pc_fhhos4_krylov fcg -pc_fhhos4_max_it 100 -pc_fhhos4_rtol 1e-10");
	EXPECT_EQ(result.ErrorCode, PETSC_SUCCESS) << LastError;
	EXPECT_LE(result.Error, 1e-8);
}

// The operator as one AIJ matrix, its rows shuffled on each process, with the index sets of its cells and faces
// (PCFHHOS4SetIS()): the iterations of the MATNEST, symmetric (2 components) or not. Every fifth face of each process
// eliminated (MatZeroRowsColumns()) and left out of the index sets: the identity there, the same solution.
TEST(PETScPC, MonolithicOperator)
{
	string options = Degrees + " -pc_fhhos4_vector_components 2";
	PETScSystem nest(0, 2, true), monolithic(0, 2, true), dirichlet(0, 2, true);
	monolithic.MakeMonolithic(0);
	PETScSystem::Result expected = nest.Solve(KSPFCG, "nestsym_", options);
	PETScSystem::Result result = monolithic.Solve(KSPFCG, "monosym_", options);
	ExpectSolves(result, monolithic.ExpectedIterations);
	EXPECT_EQ(result.Iterations, expected.Iterations);
	PETScSystem nonsymmetricNest(0.3, 1, false), nonsymmetric(0.3, 1, false);
	nonsymmetric.MakeMonolithic(0);
	expected = nonsymmetricNest.Solve(KSPFGMRES, "nestnonsym_", Degrees);
	result = nonsymmetric.Solve(KSPFGMRES, "mononsym_", Degrees);
	ExpectSolves(result, nonsymmetric.ExpectedIterations);
	EXPECT_EQ(result.Iterations, expected.Iterations);
	dirichlet.MakeMonolithic(5);
	result = dirichlet.Solve(KSPFCG, "dirichlet_", options);
	EXPECT_EQ(result.ErrorCode, PETSC_SUCCESS) << LastError;
	EXPECT_GT(result.Reason, 0);
	EXPECT_LE(result.Error, 1e-8);
}

// The discretization inferred from the blocks (ADR-0023), only -pc_fhhos4_vector_components given: the iterations with
// all the options of the discretization given, symmetric or not
TEST(PETScPC, InferredDiscretization)
{
	PETScSystem symmetric(0, 1, true), nonsymmetric(0.3, 1, false);
	PETScSystem::Result given = symmetric.Solve(KSPFCG, "givensym_", Degrees);
	PETScSystem::Result inferred = symmetric.Solve(KSPFCG, "inferredsym_", " -pc_fhhos4_vector_components 1");
	ExpectSolves(inferred, symmetric.ExpectedIterations);
	EXPECT_EQ(inferred.Iterations, given.Iterations);
	given = nonsymmetric.Solve(KSPFGMRES, "givennonsym_", Degrees);
	inferred = nonsymmetric.Solve(KSPFGMRES, "inferrednonsym_", " -pc_fhhos4_vector_components 1");
	ExpectSolves(inferred, nonsymmetric.ExpectedIterations);
	EXPECT_EQ(inferred.Iterations, given.Iterations);
}

// An ambiguity of the inference (only -pc_fhhos4_vector_components 1): blocks of 64 DoFs per cell and 9 per face (the
// system ⊗ dense blocks: 3D Q k=2 l=3, or 2D Q k=8 l=7, whose p-levels differ), refused on all the processes, the
// message naming the option to give. With it, or with the dimension and the degrees without the space (read from the
// blocks: Q), the inference resolves: the setup goes on, to the near null space, not given here.
TEST(PETScPC, Ambiguity)
{
	int n = 4;
	CartesianHybridSystem global(n, 0);
	PartitionedSystem system(global, n, PartitionShape::Blocks, World());
	std::mt19937 generator(5);
	std::uniform_real_distribution<double> coefficient(0.5, 1);
	Eigen::MatrixXd C_TT = Eigen::MatrixXd::NullaryExpr(64, 64, [&]() { return coefficient(generator); });
	Eigen::MatrixXd C_TF = Eigen::MatrixXd::NullaryExpr(64, 9, [&]() { return coefficient(generator); });
	Eigen::MatrixXd C_FF = Eigen::MatrixXd::NullaryExpr(9, 9, [&]() { return coefficient(generator); });
	SparseMatrix globals[4] = { Kron(system.Global.A_TT, C_TT), Kron(system.Global.A_TF, C_TF), Kron(system.Global.A_FT, Eigen::MatrixXd(C_TF.transpose())), Kron(system.Global.A_FF, C_FF) };
	Mat blocks[4], P = Nest(globals, *system.Cells, *system.Faces, 64, 9, blocks);
	Vec b, x;
	PetscCallAbort(PETSC_COMM_WORLD, MatCreateVecs(P, &x, &b));
	PetscCallAbort(PETSC_COMM_WORLD, VecSet(b, 1));
	// The error of the setup with these options (prefix), its message in LastError
	auto setUp = [&](const string& prefix, const string& options)
	{
		KSP ksp;
		PC pc;
		PetscCallAbort(PETSC_COMM_WORLD, KSPCreate(PETSC_COMM_WORLD, &ksp));
		PetscCallAbort(PETSC_COMM_WORLD, KSPSetOptionsPrefix(ksp, prefix.c_str()));
		PetscCallAbort(PETSC_COMM_WORLD, KSPSetOperators(ksp, P, P));
		PetscCallAbort(PETSC_COMM_WORLD, KSPSetType(ksp, KSPFCG));
		PetscCallAbort(PETSC_COMM_WORLD, KSPGetPC(ksp, &pc));
		PetscCallAbort(PETSC_COMM_WORLD, PCSetType(pc, PCFHHOS4));
		istringstream in(options);
		string token, prefixed;
		while (in >> token)
			prefixed += (token[0] == '-' ? "-" + prefix + token.substr(1) : token) + " ";
		PetscCallAbort(PETSC_COMM_WORLD, PetscOptionsInsertString(nullptr, prefixed.c_str()));
		PetscCallAbort(PETSC_COMM_WORLD, KSPSetFromOptions(ksp));
		PetscCallAbort(PETSC_COMM_WORLD, PetscPushErrorHandler(StoreError, nullptr));
		LastError.clear();
		PetscErrorCode error = KSPSolve(ksp, b, x);
		PetscCallAbort(PETSC_COMM_WORLD, PetscPopErrorHandler());
		PetscCallAbort(PETSC_COMM_WORLD, KSPDestroy(&ksp));
		return error;
	};
	EXPECT_NE(setUp("ambiguity_", "-pc_fhhos4_vector_components 1"), PETSC_SUCCESS);
	EXPECT_EQ(LastError, "fhhos4: 64 DoFs per cell and 9 per face (1 component) fit discretizations that need different solvers:\n"
	                     "  Dimension 3, Q, face degree 2, cell degree 3\n"
	                     "  Dimension 2, Q, face degree 8, cell degree 7\n"
	                     "Set -pc_fhhos4_dimension (2 or 3).");
	for (auto [prefix, options] : { pair<string, string>{ "dimension_", "-pc_fhhos4_vector_components 1 -pc_fhhos4_dimension 3" },
	                                { "degrees_", "-pc_fhhos4_dimension 3 -pc_fhhos4_face_degree 2 -pc_fhhos4_cell_degree 3" } })
	{
		EXPECT_NE(setUp(prefix, options), PETSC_SUCCESS);
		EXPECT_NE(LastError.find("near null space"), string::npos) << options << ": " << LastError;
	}
	// 2 uncoupled components, P k = l = 1 (6 DoFs per cell, 4 per face), the zeros between them not stored: the stored
	// entries give no discretization, so with the dimension and the degrees given, the space is P, as in Solver::Setup()
	auto uncoupled = [&](int rows, int cols)
	{
		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(2 * rows, 2 * cols);
		for (int component = 0; component < 2; component++)
			C.block(component * rows, component * cols, rows, cols) = Eigen::MatrixXd::NullaryExpr(rows, cols, [&]() { return coefficient(generator); });
		return C;
	};
	Eigen::MatrixXd U_TF = uncoupled(3, 2);
	SparseMatrix uncoupledGlobals[4] = { Kron(system.Global.A_TT, uncoupled(3, 3)), Kron(system.Global.A_TF, U_TF), Kron(system.Global.A_FT, Eigen::MatrixXd(U_TF.transpose())), Kron(system.Global.A_FF, uncoupled(2, 2)) };
	PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&x));
	PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&b));
	PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&P));
	for (Mat& A : blocks)
		PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&A));
	P = Nest(uncoupledGlobals, *system.Cells, *system.Faces, 6, 4, blocks);
	PetscCallAbort(PETSC_COMM_WORLD, MatCreateVecs(P, &x, &b));
	PetscCallAbort(PETSC_COMM_WORLD, VecSet(b, 1));
	EXPECT_NE(setUp("spacep_", "-pc_fhhos4_vector_components 2 -pc_fhhos4_dimension 2 -pc_fhhos4_face_degree 1 -pc_fhhos4_cell_degree 1"), PETSC_SUCCESS);
	EXPECT_NE(LastError.find("near null space"), string::npos) << LastError;
	PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&x));
	PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&b));
	PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&P));
	for (Mat& A : blocks)
		PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&A));
}

// Refused on all the processes: the degrees not given (inferred) without -pc_fhhos4_vector_components, no near null
// space, a near null space that is not the interpolation of 1 (the constant vector on 2 components), the DoFs of the
// processes not whole blocks (7 components: neither the 3600 cells nor the 7080 faces are multiples of 7), a monolithic
// operator without PCFHHOS4SetIS()
TEST(PETScPC, Errors)
{
	PETScSystem system(0, 1, true);
	EXPECT_NE(system.Solve(KSPFCG, "nodegree_", " -pc_fhhos4_dimension 2").ErrorCode, PETSC_SUCCESS);
	EXPECT_NE(LastError.find("give -pc_fhhos4_vector_components"), string::npos) << LastError;
	PETScSystem noNearNullSpace(0, 1, true, NearNullSpace::None);
	EXPECT_NE(noNearNullSpace.Solve(KSPFCG, "nonullspace_", Degrees).ErrorCode, PETSC_SUCCESS);
	PETScSystem constant(0, 2, true, NearNullSpace::Constant);
	EXPECT_NE(constant.Solve(KSPFCG, "constant_", Degrees + " -pc_fhhos4_vector_components 2").ErrorCode, PETSC_SUCCESS);
	EXPECT_NE(system.Solve(KSPFCG, "blocksize_", Degrees + " -pc_fhhos4_vector_components 7").ErrorCode, PETSC_SUCCESS);
	Mat monolithic;
	PetscCallAbort(PETSC_COMM_WORLD, MatConvert(system.P, MATAIJ, MAT_INITIAL_MATRIX, &monolithic));
	std::swap(monolithic, system.P);
	EXPECT_NE(system.Solve(KSPFCG, "monolithic_", Degrees).ErrorCode, PETSC_SUCCESS);
	EXPECT_NE(LastError.find("PCFHHOS4SetIS()"), string::npos) << LastError;
	std::swap(monolithic, system.P);
	PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&monolithic));
}

int main(int argc, char** argv)
{
	MPIEnvironment mpi(argc, argv);
	const Communicator& world = mpi.World();
	CurrentCommunicator current(world);
	::testing::InitGoogleTest(&argc, argv);
	if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS) // -dll_append: registers PCFHHOS4
		return 1;
	if (!world.IsRoot())
	{
		::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
		delete listeners.Release(listeners.default_result_printer());
		cout.rdbuf(nullptr);
	}
	int result = RUN_ALL_TESTS();
	PetscFinalize();
	return world.Max(result);
}
