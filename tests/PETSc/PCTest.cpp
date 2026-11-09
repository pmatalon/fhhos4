// fhhos4's preconditioner for PETSc (library/petsc, -pc_type fhhos4), loaded by -dll_append as by any PETSc program
// (ctest PETScPC_np1...4): the Cartesian hybrid system partitioned among the processes (PartitionedSystem), each
// process giving the rows of its cells and faces to a MATNEST of AIJ blocks (MATSEQAIJ on one process, MATMPIAIJ on
// several), solved by PETSc's Krylov methods preconditioned by fhhos4, against the library's solve of the whole system
// on one process. The results printed by the first process; a test fails if it fails on one process.
#include <sstream>
#include <gtest/gtest.h>
#include <petscksp.h>
#include <fhhos4>
#include "fhhos4_petsc.h"
#include "support/CartesianHybridSystem.h"
#include "support/Kronecker.h"
#include "support/PartitionedSystem.h"

using namespace fhhos4_tests;

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

	// The parts of a vector of the operator on the cells and faces of the process
	void SetParts(Mat P, Vec v, const Vector& v_T, const Vector& v_F)
	{
		IS rows[2];
		PetscCallAbort(PETSC_COMM_WORLD, MatNestGetISs(P, rows, nullptr));
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

	Vector Part(Mat P, Vec v, int i)
	{
		IS rows[2];
		Vec part;
		const PetscScalar* a;
		PetscInt n;
		PetscCallAbort(PETSC_COMM_WORLD, MatNestGetISs(P, rows, nullptr));
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
	// for the cells, one for the faces)
	struct PETScSystem
	{
		int c;
		Mat Blocks[4] = {}, P = nullptr;
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
			Vector cellOnes = Vector::Ones(global.A_TT.rows()), faceOnes = Vector::Ones(global.S.rows());
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
			const IndexMap* maps[4][2] = { { &cells, &cells }, { &cells, &faces }, { &faces, &cells }, { &faces, &faces } };
			const SparseMatrix* globals[4] = { &system.Global.A_TT, &system.Global.A_TF, &system.Global.A_FT, &system.Global.A_FF };
			for (int k = 0; k < 4; k++)
			{
				SparseMatrix rows = OnEachComponent(*globals[k], c).middleRows(c * maps[k][0]->Begin, c * maps[k][0]->NOwned);
				Blocks[k] = AIJ(rows, c * maps[k][1]->NOwned);
			}
			PetscCallAbort(PETSC_COMM_WORLD, MatCreateNest(PETSC_COMM_WORLD, 2, nullptr, 2, nullptr, Blocks, &P));
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
			SetParts(P, b, b_T.segment(c * cells.Begin, c * cells.NOwned), b_F.segment(c * faces.Begin, c * faces.NOwned));
			if (nearNullSpace != NearNullSpace::None)
			{
				vector<Vec> ones(nearNullSpace == NearNullSpace::Constant ? 0 : c);
				for (int k = 0; k < (int)ones.size(); k++)
				{
					// The component k of 1: the first function of each block on that component
					auto onComponent = [&](BigNumber blocks) { Vector v = Vector::Zero(c * blocks); v(Eigen::seqN(k, blocks, c)).setOnes(); return v; };
					PetscCallAbort(PETSC_COMM_WORLD, VecDuplicate(b, &ones[k]));
					SetParts(P, ones[k], onComponent(cells.NOwned), onComponent(faces.NOwned));
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
			PetscCallAbort(PETSC_COMM_WORLD, VecDestroy(&b));
			PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&P));
			for (Mat& A : Blocks)
				PetscCallAbort(PETSC_COMM_WORLD, MatDestroy(&A));
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
			PetscCallAbort(PETSC_COMM_WORLD, PetscPushErrorHandler(PetscReturnErrorHandler, nullptr));
			result.ErrorCode = KSPSolve(ksp, b, x);
			PetscCallAbort(PETSC_COMM_WORLD, PetscPopErrorHandler());
			if (result.ErrorCode == PETSC_SUCCESS)
			{
				PetscCallAbort(PETSC_COMM_WORLD, KSPGetIterationNumber(ksp, &result.Iterations));
				PetscCallAbort(PETSC_COMM_WORLD, KSPGetConvergedReason(ksp, &result.Reason));
				Vector allX(ExpectedX_T.rows() + ExpectedX_F.rows()), expected(allX.rows());
				allX << Gathered(Part(P, x, 0)), Gathered(Part(P, x, 1));
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
		EXPECT_EQ(result.ErrorCode, PETSC_SUCCESS);
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
	EXPECT_EQ(result.ErrorCode, PETSC_SUCCESS);
	EXPECT_LE(result.Error, 1e-8);
}

// Refused on all the processes: the degrees not given, no near null space, a near null space that is not the
// interpolation of 1 (the constant vector on 2 components), the DoFs of the processes not whole blocks (7 components:
// neither the 3600 cells nor the 7080 faces are multiples of 7), a monolithic operator
TEST(PETScPC, Errors)
{
	PETScSystem system(0, 1, true);
	EXPECT_NE(system.Solve(KSPFCG, "nodegree_", " -pc_fhhos4_dimension 2").ErrorCode, PETSC_SUCCESS);
	PETScSystem noNearNullSpace(0, 1, true, NearNullSpace::None);
	EXPECT_NE(noNearNullSpace.Solve(KSPFCG, "nonullspace_", Degrees).ErrorCode, PETSC_SUCCESS);
	PETScSystem constant(0, 2, true, NearNullSpace::Constant);
	EXPECT_NE(constant.Solve(KSPFCG, "constant_", Degrees + " -pc_fhhos4_vector_components 2").ErrorCode, PETSC_SUCCESS);
	EXPECT_NE(system.Solve(KSPFCG, "blocksize_", Degrees + " -pc_fhhos4_vector_components 7").ErrorCode, PETSC_SUCCESS);
	Mat monolithic;
	PetscCallAbort(PETSC_COMM_WORLD, MatConvert(system.P, MATAIJ, MAT_INITIAL_MATRIX, &monolithic));
	std::swap(monolithic, system.P);
	EXPECT_NE(system.Solve(KSPFCG, "monolithic_", Degrees).ErrorCode, PETSC_SUCCESS);
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
