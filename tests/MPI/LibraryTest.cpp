#include <random>
#include <gtest/gtest.h>
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <fhhos4>
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

	fhhos4::Solver CartesianSolver(MPI_Comm comm)
	{
		fhhos4::Solver solver;
		solver.Comm = comm;
		solver.Dimension = 2;
		solver.FaceDegree = 0;
		solver.CellDegree = 0;
		solver.Tolerance = 1e-10;
		return solver;
	}

	// The vectors of all the processes, in rank order
	Vector Gathered(const Vector& local)
	{
		vector<double> values(local.data(), local.data() + local.rows());
		vector<double> all = World().AllGatherV(values);
		return Eigen::Map<Vector>(all.data(), all.size());
	}
}

// The entities of the columns of the rows not owned in owners (the ghosts), in reverse order of their global numbers
// (any order)
vector<BigNumber> ReversedGhosts(std::initializer_list<const SparseMatrix*> rows, const IndexMap& owners)
{
	vector<BigNumber> ghosts;
	for (const SparseMatrix* M : rows)
		for (BigNumber p = 0; p < (BigNumber)M->nonZeros(); p++)
			if (!owners.IsOwned(M->innerIndexPtr()[p]))
				ghosts.push_back(M->innerIndexPtr()[p]);
	sort(ghosts.begin(), ghosts.end());
	ghosts.erase(unique(ghosts.begin(), ghosts.end()), ghosts.end());
	reverse(ghosts.begin(), ghosts.end());
	return ghosts;
}

// The rows with their columns numbered locally: the entities owned in owners, then the ghosts in their order
SparseMatrix LocalColumns(const SparseMatrix& rows, const IndexMap& owners, const vector<BigNumber>& ghosts)
{
	auto local = [&](BigNumber e) -> BigNumber
	{
		if (owners.IsOwned(e))
			return e - owners.Begin;
		auto it = find(ghosts.begin(), ghosts.end(), e);
		return it == ghosts.end() ? -1 : owners.NOwned + (it - ghosts.begin());
	};
	return DistributedMatrix::RenumberColumns(rows, local, owners.NOwned + ghosts.size(), 1);
}

// Rows with global columns in two parts, as PETSc and hypre store them: Diag, the columns [begin, end) numbered from 0,
// and Offd, the others, compressed, whose global numbers are OffdColumns (increasing)
struct TwoParts
{
	fhhos4::SparseMatrix Diag, Offd;
	vector<int64_t> OffdColumns;
};

TwoParts Split(const fhhos4::SparseMatrix& rows, int64_t begin, int64_t end)
{
	TwoParts parts;
	for (BigNumber p = 0; p < (BigNumber)rows.nonZeros(); p++)
		if (rows.innerIndexPtr()[p] < begin || rows.innerIndexPtr()[p] >= end)
			parts.OffdColumns.push_back(rows.innerIndexPtr()[p]);
	sort(parts.OffdColumns.begin(), parts.OffdColumns.end());
	parts.OffdColumns.erase(unique(parts.OffdColumns.begin(), parts.OffdColumns.end()), parts.OffdColumns.end());
	vector<Eigen::Triplet<double>> diag, offd;
	for (BigNumber i = 0; i < (BigNumber)rows.rows(); i++)
	{
		for (fhhos4::SparseMatrix::InnerIterator it(rows, i); it; ++it)
		{
			if (it.col() >= begin && it.col() < end)
				diag.emplace_back(i, it.col() - begin, it.value());
			else
				offd.emplace_back(i, lower_bound(parts.OffdColumns.begin(), parts.OffdColumns.end(), it.col()) - parts.OffdColumns.begin(), it.value());
		}
	}
	parts.Diag.resize(rows.rows(), end - begin);
	parts.Diag.setFromTriplets(diag.begin(), diag.end());
	parts.Offd.resize(rows.rows(), parts.OffdColumns.size());
	parts.Offd.setFromTriplets(offd.begin(), offd.end());
	parts.Diag.makeCompressed();
	parts.Offd.makeCompressed();
	return parts;
}

// The uncondensed system distributed by rows, through the interface of the library: the ghost faces and cells numbered
// in reverse order of their global numbers (any order). The solution of the system on one process. Non-symmetric if
// asymmetry != 0. singular: the pure Neumann system (SingularMatrix; the right-hand side incompatible: projected),
// the solutions compared up to the kernel (the constants, of each component). components: the system on each of them
// (the blocks ⊗ I_n, VectorComponents = n; the ghosts are the same faces and cells). mixedPrecision, blockCSR:
// MixedPrecision and BlockCSRStorage of both solvers (the defaults: true, false). mergeOffDiag: the blocks of each
// process in two parts (Split(), as PETSc stores them), merged by MergeOffDiag() (the ghosts
// by increasing global number), checked against the test's. shape: the partition. inferred: the distributed solver with
// only VectorComponents set (the discretization inferred from the blocks of the processes), checked after the setup.
// Returns the iterations on the processes; setupLog: the printouts of their setup (Verbosity 2, on the root process)
// if not null.
int ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation prolongation, double asymmetry = 0, bool singular = false, int components = 1, bool mixedPrecision = true, string* setupLog = nullptr, bool blockCSR = false, bool mergeOffDiag = false, PartitionShape shape = PartitionShape::Blocks, bool inferred = false)
{
	int n = 60;
	CartesianHybridSystem global(n, asymmetry, singular);
	PartitionedSystem system(global, n, shape, World());
	const IndexMap& faces = *system.Faces;
	const IndexMap& cells = *system.Cells;
	int c = components;
	Vector b_T = OnEachComponent(Vector::LinSpaced(global.A_TT.rows(), 1, 2), c), b_F = Vector::Zero(c * global.S.rows());

	// One process: the whole system
	Vector expectedX_T(b_T.rows()), expectedX_F = Vector::Zero(b_F.rows());
	Vector allX = Vector::LinSpaced(b_F.rows(), -1, 2), expectedAx(b_F.rows()); // ApplyCondensedMatrix()
	int expectedIterations;
	{
		fhhos4::Solver solver = CartesianSolver(MPI_COMM_SELF);
		solver.Prolongation = prolongation;
		solver.SingularMatrix = singular;
		solver.VectorComponents = c;
		solver.MixedPrecision = mixedPrecision;
		solver.BlockCSRStorage = blockCSR;
		Eigen::MatrixXd cellOnes = Constants(global.A_TT.rows(), c), faceOnes = Constants(global.S.rows(), c);
		fhhos4::SparseMatrix A_TT = OnEachComponent(system.Global.A_TT, c), A_TF = OnEachComponent(system.Global.A_TF, c);
		fhhos4::SparseMatrix A_FT = OnEachComponent(system.Global.A_FT, c), A_FF = OnEachComponent(system.Global.A_FF, c);
		if (asymmetry == 0)
			solver.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes);
		else
			solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellOnes, faceOnes);
		expectedIterations = solver.Solve(b_T, b_F, expectedX_T, expectedX_F).Iterations;
		solver.ApplyCondensedMatrix(allX, expectedAx);
	}

	// The faces of the rows of the process (of A_TF and A_FF): its own, then the others; the cells of A_FT: the same
	SparseMatrix ownA_TF = system.Global.A_TF.middleRows(cells.Begin, cells.NOwned);
	SparseMatrix ownA_FF = system.Global.A_FF.middleRows(faces.Begin, faces.NOwned);
	SparseMatrix ownA_FT = system.Global.A_FT.middleRows(faces.Begin, faces.NOwned);
	vector<BigNumber> ghostFaces = ReversedGhosts({ &ownA_TF, &ownA_FF }, faces), ghostCells = ReversedGhosts({ &ownA_FT }, cells);
	fhhos4::SparseMatrix A_TT = OnEachComponent(system.A_TT, c), A_TF = OnEachComponent(LocalColumns(ownA_TF, faces, ghostFaces), c);
	fhhos4::SparseMatrix A_FF = OnEachComponent(LocalColumns(ownA_FF, faces, ghostFaces), c), A_FT = OnEachComponent(LocalColumns(ownA_FT, cells, ghostCells), c);
	vector<int64_t> faceGhostGlobalNumbers(ghostFaces.begin(), ghostFaces.end()), cellGhostGlobalNumbers(ghostCells.begin(), ghostCells.end());
	if (mergeOffDiag)
	{
		// The rows with global DoF columns (the DoF k of the face f: c f + k)
		int64_t faceBegin = c * faces.Begin, faceEnd = c * (faces.Begin + faces.NOwned), cellBegin = c * cells.Begin, cellEnd = c * (cells.Begin + cells.NOwned);
		TwoParts tf = Split(OnEachComponent(ownA_TF, c), faceBegin, faceEnd), ff = Split(OnEachComponent(ownA_FF, c), faceBegin, faceEnd);
		TwoParts ft = Split(OnEachComponent(ownA_FT, c), cellBegin, cellEnd);
		fhhos4::DistributedBlocks blocks = asymmetry == 0
			? fhhos4::MergeOffDiag(tf.Diag, tf.Offd, tf.OffdColumns.data(), ff.Diag, ff.Offd, ff.OffdColumns.data(), c)
			: fhhos4::MergeOffDiag(tf.Diag, tf.Offd, tf.OffdColumns.data(), ft.Diag, ft.Offd, ft.OffdColumns.data(), ff.Diag, ff.Offd, ff.OffdColumns.data(), c, c);
		sort(ghostFaces.begin(), ghostFaces.end());
		sort(ghostCells.begin(), ghostCells.end());
		faceGhostGlobalNumbers.assign(ghostFaces.begin(), ghostFaces.end());
		cellGhostGlobalNumbers = asymmetry == 0 ? vector<int64_t>() : vector<int64_t>(ghostCells.begin(), ghostCells.end());
		EXPECT_EQ(blocks.FaceGhostGlobalNumbers, faceGhostGlobalNumbers);
		EXPECT_EQ(blocks.CellGhostGlobalNumbers, cellGhostGlobalNumbers);
		EXPECT_EQ(SparseMatrix(blocks.A_TF - OnEachComponent(LocalColumns(ownA_TF, faces, ghostFaces), c)).norm(), 0);
		EXPECT_EQ(SparseMatrix(blocks.A_FF - OnEachComponent(LocalColumns(ownA_FF, faces, ghostFaces), c)).norm(), 0);
		if (asymmetry == 0)
			EXPECT_EQ(blocks.A_FT.rows(), 0);
		else
			EXPECT_EQ(SparseMatrix(blocks.A_FT - OnEachComponent(LocalColumns(ownA_FT, cells, ghostCells), c)).norm(), 0);
		A_TF = blocks.A_TF;
		A_FF = blocks.A_FF;
		A_FT = blocks.A_FT;
	}

	fhhos4::Solver solver = CartesianSolver(World().Comm());
	solver.Prolongation = prolongation;
	solver.SingularMatrix = singular;
	solver.VectorComponents = c;
	solver.MixedPrecision = mixedPrecision;
	solver.BlockCSRStorage = blockCSR;
	Eigen::MatrixXd cellOnes = Constants(cells.NOwned, c), faceOnes = Constants(faces.NOwned, c);
	if (setupLog)
	{
		solver.Verbosity = 2;
		testing::internal::CaptureStdout();
	}
	if (inferred)
	{
		solver.Dimension = -1;
		solver.FaceDegree = -1;
		solver.CellDegree = -1;
	}
	if (asymmetry == 0)
		solver.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes, faceGhostGlobalNumbers);
	else
		solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellOnes, faceOnes, faceGhostGlobalNumbers, cellGhostGlobalNumbers);
	EXPECT_EQ(solver.FaceBlockSize(), c);
	EXPECT_EQ(solver.CellBlockSize(), c);
	if (setupLog)
	{
		*setupLog = testing::internal::GetCapturedStdout();
		solver.Verbosity = 0;
	}
	Vector x_T(c * cells.NOwned), x_F = Vector::Zero(c * faces.NOwned);
	fhhos4::Result result = solver.Solve(b_T.segment(c * cells.Begin, c * cells.NOwned), b_F.segment(c * faces.Begin, c * faces.NOwned), x_T, x_F);

	EXPECT_TRUE(result.Converged);
	// singular: with HarmonicExtension, 14 iterations on one process (and on 3), 18-19 on 2 and 4, as the regular system
	EXPECT_LE(result.Iterations, expectedIterations + (singular ? 6 : 4));
	if (World().Size() == 1)
		EXPECT_EQ(result.Iterations, expectedIterations);
	Vector allX_T = Gathered(x_T), allX_F = Gathered(x_F);
	if (singular)
	{
		// The difference without its component along the kernel (the interpolation of 1: ones, on each component: the
		// rows of the vectors reshaped by component)
		Vector shift = ((allX_T - expectedX_T).reshaped(c, allX_T.rows() / c).rowwise().sum() + (allX_F - expectedX_F).reshaped(c, allX_F.rows() / c).rowwise().sum()) / double((allX_T.rows() + allX_F.rows()) / c);
		allX_T.reshaped(c, allX_T.rows() / c).colwise() -= shift;
		allX_F.reshaped(c, allX_F.rows() / c).colwise() -= shift;
	}
	EXPECT_LE((allX_F - expectedX_F).norm(), 1e-8 * expectedX_F.norm());
	EXPECT_LE((allX_T - expectedX_T).norm(), 1e-8 * expectedX_T.norm());
	// The product by the condensed matrix of the processes: that of the whole system
	Vector Ax(c * faces.NOwned);
	solver.ApplyCondensedMatrix(allX.segment(c * faces.Begin, c * faces.NOwned), Ax);
	EXPECT_LE((Gathered(Ax) - expectedAx).norm(), 1e-13 * expectedAx.norm());
	if (World().IsRoot())
		cout << "\t" << World().Size() << " processes: " << result.Iterations << " iterations (" << expectedIterations << " on one)" << endl;
	return result.Iterations;
}

TEST(DistributedLibrary, UncondensedSystem)
{
	ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::Reconstruction);
}

// Prolongation::HarmonicExtension: the local solves of the removed faces in the aggregates of each process
TEST(DistributedLibrary, UncondensedSystemHarmonicExtension)
{
	ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::HarmonicExtension);
}

// SetupNonSym(): A_FT with the ghost cells of its rows; 0.9: convection-dominated (the downwind order of each
// process)
TEST(DistributedLibrary, NonsymmetricUncondensedSystem)
{
	for (double asymmetry : { 0.3, 0.9 })
	{
		ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::Reconstruction, asymmetry);
		ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::HarmonicExtension, asymmetry);
	}
}

// A vector-valued unknown (VectorComponents): the blocks of each process ⊗ I_n, the ghosts unchanged, symmetric or not,
// with both prolongations; singular: a pinned DoF per component on the coarsest level
TEST(DistributedLibrary, VectorComponents)
{
	ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::HarmonicExtension, 0, false, 2);
	ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::HarmonicExtension, 0.3, false, 3);
	ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::Reconstruction, 0, false, 2);
	ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::HarmonicExtension, 0, true, 2);
}

// SingularMatrix: the pinned DoF of the coarsest level on one process, its column among the ghosts of the others
TEST(DistributedLibrary, SingularMatrix)
{
	ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::Reconstruction, 0, true);
	ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::HarmonicExtension, 0, true);
}

// Mixed precision (the default: the levels in single precision, in block CSR with 2 components) and in double: the same
// iterations within 1, as on one process (UAMGMixedPrecisionProgram); the levels in single precision on the processes.
// Symmetric (the Chebyshev smoothers) and not (asymmetry 0.3: the l1 hybrid block Gauss-Seidel).
TEST(DistributedLibrary, MixedPrecision)
{
	for (double asymmetry : { 0.0, 0.3 })
	{
		for (int components : { 1, 2 })
		{
			string setupLog;
			int mixed = ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::Reconstruction, asymmetry, false, components, true, &setupLog);
			int inDouble = ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::Reconstruction, asymmetry, false, components, false);
			EXPECT_LE(abs(mixed - inDouble), 1) << components << " components, asymmetry " << asymmetry;
			if (World().IsRoot())
			{
				EXPECT_EQ(setupLog.find("Levels in double precision"), string::npos) << setupLog;
				EXPECT_NE(setupLog.find("single precision)"), string::npos) << setupLog;
				if (components == 2)
					EXPECT_NE(setupLog.find("single precision, in block CSR (blocks of 2))"), string::npos) << setupLog;
			}
		}
	}
}

// The block CSR storage in double (BlockCSRStorage; 2 components: blocks of 2 on all the levels but the coarsest) and
// CSR: the same iterations on 1-4 processes, the levels in block CSR on the processes. Symmetric (the Chebyshev
// smoothers) and not (the l1 hybrid block Gauss-Seidel).
TEST(DistributedLibrary, BlockCSR)
{
	for (double asymmetry : { 0.0, 0.3 })
	{
		string setupLog;
		int blocks = ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::Reconstruction, asymmetry, false, 2, false, &setupLog, true);
		int csr = ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::Reconstruction, asymmetry, false, 2, false);
		EXPECT_EQ(blocks, csr) << "asymmetry " << asymmetry;
		if (World().IsRoot())
			EXPECT_NE(setupLog.find("(and in block CSR (blocks of 2))"), string::npos) << setupLog;
	}
}

// MergeOffDiag(): the blocks of each process in two parts, as PETSc and hypre store them (the
// columns of its own faces or cells, and the others with their global numbers), merged; on 1 and 2 components
TEST(DistributedLibrary, MergeOffDiag)
{
	for (double asymmetry : { 0.0, 0.3 })
		for (int components : { 1, 2 })
			ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::HarmonicExtension, asymmetry, false, components, true, nullptr, false, true);
}

// A process that owns no cell and no face (NGSolve's rank 0): its levels empty, in block CSR (2 components, mixed
// precision) too; symmetric (the Chebyshev smoothers) and not (the l1 hybrid block Gauss-Seidel)
TEST(DistributedLibrary, EmptyProcess)
{
	if (World().Size() == 1)
		GTEST_SKIP() << "several processes";
	for (double asymmetry : { 0.0, 0.3 })
		ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::HarmonicExtension, asymmetry, false, 2, true, nullptr, false, false, PartitionShape::StripsWithFirstProcessEmpty);
}

// The discretization inferred from the blocks of the processes (only VectorComponents set; ADR-0023), on blocks of
// processes and with a process that owns no cell (no opinion: it receives the result), symmetric or not: the
// iterations with all the parameters set
TEST(DistributedLibrary, InferredDiscretization)
{
	for (double asymmetry : { 0.0, 0.3 })
	{
		for (PartitionShape shape : { PartitionShape::Blocks, PartitionShape::StripsWithFirstProcessEmpty })
		{
			if (shape == PartitionShape::StripsWithFirstProcessEmpty && World().Size() == 1)
				continue;
			int set = ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::HarmonicExtension, asymmetry, false, 1, true, nullptr, false, false, shape);
			int inferred = ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation::HarmonicExtension, asymmetry, false, 1, true, nullptr, false, false, shape, true);
			EXPECT_EQ(inferred, set) << "asymmetry " << asymmetry;
		}
	}
}

// More near null space functions than the constants (the rigid body modes of elasticity; here, as in
// UAMGLibrary.NearNullSpaceBeyondTheConstants, a third function on 2 components): the Cartesian system ⊗ dense blocks
// (2D, k = l = 1: 6 DoFs per cell, 4 per face) whose near kernel the 3 functions span, the functions of the cells and
// faces of each process, on blocks of processes and with a process that owns nothing. With the level of the functions
// as the coarsest one (CoarseMatrixMaxRows), the multigrid has no aggregation: that of one process, and its iterations
// within 1 (the operator of that level reads the functions of the ghost faces, received from their owners: wrong ones
// cost 83 iterations instead of 5); with the h-levels under it, each process aggregates its cells. Singular
// (SingularMatrix: the boundary faces in the system, whose kernel the 3 functions then are, with their exact cells
// -A_TT^-1 A_TF, mixed so that they are not orthogonal): the same, and the finest level as the coarsest, pinned on the
// mixed functions. Refused on all the processes: dependent functions on a cell of one process, another number of
// functions on one process.
TEST(DistributedLibrary, NearNullSpaceBeyondTheConstants)
{
	int n = 40;
	CartesianHybridSystem global(n, 0), neumann(n, 0, true);
	// The dense blocks of DenseBlocks() (UAMGLibraryTest.cpp): C_TT symmetric positive definite, G = C_TF^T C_TT^-1 C_TF
	std::mt19937 generator(6);
	std::uniform_real_distribution<double> coefficient(0.5, 1);
	Eigen::MatrixXd B = Eigen::MatrixXd::NullaryExpr(6, 6, [&]() { return coefficient(generator); });
	Eigen::MatrixXd C_TT = B * B.transpose() + 6 * Eigen::MatrixXd::Identity(6, 6);
	Eigen::MatrixXd C_TF = Eigen::MatrixXd::NullaryExpr(6, 4, [&]() { return coefficient(generator); });
	Eigen::MatrixXd G = C_TF.transpose() * C_TT.llt().solve(C_TF);
	// Added to the block of each face: the DoF 3 out of the near kernel, which the 3 functions then span (S ⊗ G alone
	// has all the 1 ⊗ v)
	Eigen::MatrixXd H = Eigen::MatrixXd::Zero(4, 4);
	H(3, 3) = 4 * G.diagonal().maxCoeff();
	// The near kernel on the faces scaled by d (the faces of the system scaled by d_f, different from one face to the
	// next: the functions of the ghosts must be those of their faces): 1/d_f on the DoFs 0, 1 (the constants of the 2
	// components) and 2 of each face; on the cells, those DoFs
	auto functions = [](const Vector& scale, int blockSize)
	{
		Eigen::MatrixXd F = Eigen::MatrixXd::Zero(scale.rows() * blockSize, 3);
		for (Eigen::Index b = 0; b < scale.rows(); b++)
			F.block(b * blockSize, 0, 3, 3) = Eigen::MatrixXd::Identity(3, 3) / scale[b];
		return F;
	};
	struct Solution
	{
		int Iterations;
		Vector X_T, X_F;
	};
	auto solve = [](MPI_Comm comm, int coarseMatrixMaxRows, const fhhos4::SparseMatrix& A_TT, const fhhos4::SparseMatrix& A_TF, const fhhos4::SparseMatrix& A_FF,
	                const Eigen::MatrixXd& cellFunctions, const Eigen::MatrixXd& faceFunctions, const vector<int64_t>& ghosts, const Vector& b_T, const Vector& b_F, bool singular = false)
	{
		fhhos4::Solver solver = CartesianSolver(comm);
		solver.SingularMatrix = singular;
		solver.FaceDegree = 1;
		solver.CellDegree = 1;
		solver.VectorComponents = 2;
		if (coarseMatrixMaxRows > 0)
			solver.CoarseMatrixMaxRows = coarseMatrixMaxRows;
		solver.SetupSym(A_TT, A_TF, A_FF, cellFunctions, faceFunctions, ghosts);
		Solution solution{ 0, Vector(b_T.rows()), Vector::Zero(b_F.rows()) };
		fhhos4::Result result = solver.Solve(b_T, b_F, solution.X_T, solution.X_F);
		EXPECT_TRUE(result.Converged);
		solution.Iterations = result.Iterations;
		return solution;
	};

	for (bool singular : { false, true })
	for (PartitionShape shape : { PartitionShape::Blocks, PartitionShape::StripsWithFirstProcessEmpty })
	{
		if (shape == PartitionShape::StripsWithFirstProcessEmpty && World().Size() == 1)
			continue;
		PartitionedSystem system(singular ? neumann : global, n, shape, World());
		const IndexMap& faces = *system.Faces;
		const IndexMap& cells = *system.Cells;
		Vector d = Vector::LinSpaced(faces.GlobalSize, 0, faces.GlobalSize - 1).unaryExpr([](double f) { return 1 + fmod(f, 5) / 4; });
		Eigen::MatrixXd cellFunctions = functions(Vector::Ones(cells.GlobalSize), 6);
		Eigen::MatrixXd faceFunctions = functions(d, 4);
		Vector b_T = Vector::LinSpaced(6 * cells.GlobalSize, 1, 2);
		Vector b_F = Vector::Zero(4 * faces.GlobalSize);
		int functionsLevelRows = 3 * (int)faces.GlobalSize; // the coarsest: the level of the functions
		SparseMatrix scaledA_TF = system.Global.A_TF * d.asDiagonal();
		SparseMatrix scaledA_FF = d.asDiagonal() * system.Global.A_FF * d.asDiagonal();
		SparseMatrix identity(faces.GlobalSize, faces.GlobalSize);
		identity.setIdentity();
		if (singular)
		{
			Eigen::Matrix3d mix;
			mix << 1, 0, 0.5, 0.3, 1, -0.2, 0.7, 0.4, 1;
			faceFunctions *= mix;
			cellFunctions = -Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>(Eigen::SparseMatrix<double>(Kron(system.Global.A_TT, C_TT))).solve(Eigen::MatrixXd(Kron(scaledA_TF, C_TF) * faceFunctions));
		}

		SparseMatrix ownA_TF = scaledA_TF.middleRows(cells.Begin, cells.NOwned);
		SparseMatrix ownA_FF = scaledA_FF.middleRows(faces.Begin, faces.NOwned);
		SparseMatrix ownIdentity = identity.middleRows(faces.Begin, faces.NOwned);
		vector<BigNumber> ghostFaces = ReversedGhosts({ &ownA_TF, &ownA_FF }, faces);
		fhhos4::SparseMatrix A_TT = Kron(system.A_TT, C_TT);
		fhhos4::SparseMatrix A_TF = Kron(LocalColumns(ownA_TF, faces, ghostFaces), C_TF);
		fhhos4::SparseMatrix A_FF = Kron(LocalColumns(ownA_FF, faces, ghostFaces), G) + Kron(LocalColumns(ownIdentity, faces, ghostFaces), H);
		Eigen::MatrixXd ownCellFunctions = cellFunctions.middleRows(6 * cells.Begin, 6 * cells.NOwned);
		Eigen::MatrixXd ownFaceFunctions = faceFunctions.middleRows(4 * faces.Begin, 4 * faces.NOwned);
		vector<int64_t> ghosts(ghostFaces.begin(), ghostFaces.end());
		Vector ownB_T = b_T.segment(6 * cells.Begin, 6 * cells.NOwned);
		Vector ownB_F = b_F.segment(4 * faces.Begin, 4 * faces.NOwned);

		vector<int> coarseMatrixMaxRowsCases = { functionsLevelRows, 0 };
		if (singular)
			coarseMatrixMaxRowsCases.push_back(4 * (int)faces.GlobalSize); // the coarsest: the finest
		for (int coarseMatrixMaxRows : coarseMatrixMaxRowsCases)
		{
			Solution expected = solve(MPI_COMM_SELF, coarseMatrixMaxRows, Kron(system.Global.A_TT, C_TT), Kron(scaledA_TF, C_TF), Kron(scaledA_FF, G) + Kron(identity, H), cellFunctions, faceFunctions, {}, b_T, b_F, singular);
			Solution distributed = solve(World().Comm(), coarseMatrixMaxRows, A_TT, A_TF, A_FF, ownCellFunctions, ownFaceFunctions, ghosts, ownB_T, ownB_F, singular);
			string where = to_string(World().Size()) + " processes, shape " + to_string((int)shape) + ", CoarseMatrixMaxRows " + to_string(coarseMatrixMaxRows) + (singular ? ", singular" : "");
			Vector allX_T = Gathered(distributed.X_T);
			Vector allX_F = Gathered(distributed.X_F);
			EXPECT_LE((allX_F - expected.X_F).norm(), 1e-8 * expected.X_F.norm()) << where;
			EXPECT_LE((allX_T - expected.X_T).norm(), 1e-8 * expected.X_T.norm()) << where;
			if (coarseMatrixMaxRows > 0)
				EXPECT_LE(abs(distributed.Iterations - expected.Iterations), 1) << where;
			else
				EXPECT_LE(distributed.Iterations, expected.Iterations + 4) << where;
			if (World().IsRoot())
				cout << "\t" << where << ": " << distributed.Iterations << " iterations (" << expected.Iterations << " on one)" << endl;
		}

		if (singular || shape != PartitionShape::Blocks || World().Size() == 1)
			continue;
		bool last = World().Rank() == World().Size() - 1;
		Eigen::MatrixXd dependent = ownCellFunctions;
		if (last)
			dependent.block(0, 2, 6, 1) = dependent.block(0, 0, 6, 1);
		EXPECT_THROW(solve(World().Comm(), 0, A_TT, A_TF, A_FF, dependent, ownFaceFunctions, ghosts, ownB_T, ownB_F), fhhos4::Error);
		Eigen::Index nFunctions = last ? 2 : 3; // the last process: the constants only
		EXPECT_THROW(solve(World().Comm(), 0, A_TT, A_TF, A_FF, ownCellFunctions.leftCols(nFunctions), ownFaceFunctions.leftCols(nFunctions), ghosts, ownB_T, ownB_F), fhhos4::Error);
	}
}

// The solution of the uncondensed system (singular: the pure Neumann one, SingularMatrix) distributed by rows over
// blocks of processes (the ghosts in reverse order), gathered on all the processes, solved by the solver of the
// processes configured by set()
template <class Set>
Vector DistributedSolution(double asymmetry, bool singular, Set set)
{
	int n = 60;
	CartesianHybridSystem global(n, asymmetry, singular);
	PartitionedSystem system(global, n, PartitionShape::Blocks, World());
	const IndexMap& faces = *system.Faces;
	const IndexMap& cells = *system.Cells;
	SparseMatrix ownA_TF = system.Global.A_TF.middleRows(cells.Begin, cells.NOwned);
	SparseMatrix ownA_FF = system.Global.A_FF.middleRows(faces.Begin, faces.NOwned);
	SparseMatrix ownA_FT = system.Global.A_FT.middleRows(faces.Begin, faces.NOwned);
	vector<BigNumber> ghostFaces = ReversedGhosts({ &ownA_TF, &ownA_FF }, faces), ghostCells = ReversedGhosts({ &ownA_FT }, cells);
	fhhos4::SparseMatrix A_TT = system.A_TT, A_TF = LocalColumns(ownA_TF, faces, ghostFaces);
	fhhos4::SparseMatrix A_FF = LocalColumns(ownA_FF, faces, ghostFaces), A_FT = LocalColumns(ownA_FT, cells, ghostCells);
	vector<int64_t> faceGhosts(ghostFaces.begin(), ghostFaces.end()), cellGhosts(ghostCells.begin(), ghostCells.end());
	Vector cellOnes = Vector::Ones(cells.NOwned), faceOnes = Vector::Ones(faces.NOwned);
	fhhos4::Solver solver = CartesianSolver(World().Comm());
	solver.SingularMatrix = singular;
	set(solver);
	if (asymmetry == 0)
		solver.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes, faceGhosts);
	else
		solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellOnes, faceOnes, faceGhosts, cellGhosts);
	Vector b_T = Vector::LinSpaced(global.A_TT.rows(), 1, 2).segment(cells.Begin, cells.NOwned), b_F = Vector::Zero(faces.NOwned);
	Vector x_T(cells.NOwned), x_F = Vector::Zero(faces.NOwned);
	EXPECT_TRUE(solver.Solve(b_T, b_F, x_T, x_F).Converged);
	Vector x(global.A_TT.rows() + global.S.rows());
	x << Gathered(x_T), Gathered(x_F);
	return x;
}

// The reproducible modes on several processes (ADR-0024): repeated setups and solves give the same solution bit for
// bit (the coarse matrix gathered on the process 0: MUMPS's solutions vary from run to run), symmetric, singular (its
// matrix pinned before it is gathered) or not; the default smoother of OnAnyNThreads for a non-symmetric matrix, the
// l1 block Gauss-Seidel, is the l1 hybrid one bit for bit on one thread per process; CoarseSolver::MUMPS refused
TEST(DistributedLibrary, ReproducibleModes)
{
	for (fhhos4::ReproducibilityMode mode : { fhhos4::ReproducibilityMode::OnFixedNThreads, fhhos4::ReproducibilityMode::OnAnyNThreads })
	{
		for (auto [asymmetry, singular] : { pair{ 0.0, false }, { 0.3, false }, { 0.0, true } })
		{
			auto set = [mode](fhhos4::Solver& s) { s.Reproducibility = mode; };
			Vector first = DistributedSolution(asymmetry, singular, set), second = DistributedSolution(asymmetry, singular, set);
			EXPECT_TRUE((first.array() == second.array()).all()) << "mode " << (int)mode << ", asymmetry " << asymmetry << ", singular " << singular;
		}
	}

	auto oneThread = [](fhhos4::Smoother smoother) { return [smoother](fhhos4::Solver& s) { s.Threads = 1; s.PreSmoother = smoother; }; };
	Vector l1 = DistributedSolution(0.3, false, oneThread(fhhos4::Smoother::L1BlockGaussSeidel));
	Vector hybrid = DistributedSolution(0.3, false, oneThread(fhhos4::Smoother::L1HybridBlockGaussSeidel));
	EXPECT_TRUE((l1.array() == hybrid.array()).all());

	if (World().Size() == 1)
		return; // MUMPS: on one process, reproducible
	for (fhhos4::ReproducibilityMode mode : { fhhos4::ReproducibilityMode::OnFixedNThreads, fhhos4::ReproducibilityMode::OnAnyNThreads })
	{
		try
		{
			DistributedSolution(0, false, [mode](fhhos4::Solver& s) { s.Reproducibility = mode; s.CoarseSolver = fhhos4::CoarseSolver::MUMPS; });
			ADD_FAILURE() << "MUMPS accepted on several processes in a reproducible mode";
		}
		catch (const fhhos4::Error& e)
		{
			EXPECT_NE(string(e.what()).find("ReproducibilityMode::On"), string::npos) << e.what();
			EXPECT_NE(string(e.what()).find("leave CoarseSolver Default"), string::npos) << e.what();
		}
	}
}

// MinRowsPerProcess on P processes: the coarsening stops as with CoarseMatrixMaxRows at MinRowsPerProcess * P rows (on
// 2 and 4 processes, a level earlier than at 1000); at the fine level with the number of faces
TEST(DistributedLibrary, MinRowsPerProcess)
{
	if (World().Size() == 1)
		GTEST_SKIP() << "several processes";
	int n = 60;
	CartesianHybridSystem global(n, 0);
	PartitionedSystem system(global, n, PartitionShape::Blocks, World());
	fhhos4::SparseMatrix A_TT = system.A_TT, A_TF = system.A_TF.LocalColumns();
	fhhos4::SparseMatrix A_FF = DistributedMatrix::RenumberColumns(SparseMatrix(system.Global.A_FF.middleRows(system.Faces->Begin, system.Faces->NOwned)), [&](BigNumber f) { return system.A_TF.ColMap->Local(f); }, system.A_TF.ColMap->NLocal(), 1);
	vector<int64_t> ghosts(system.A_TF.ColMap->Ghosts.begin(), system.A_TF.ColMap->Ghosts.end());
	Vector cellOnes = Vector::Ones(system.Cells->NOwned), faceOnes = Vector::Ones(system.Faces->NOwned);
	auto levels = [&](int coarseMatrixMaxRows, int minRowsPerProcess)
	{
		fhhos4::Solver solver = CartesianSolver(World().Comm());
		solver.CoarseMatrixMaxRows = coarseMatrixMaxRows;
		solver.MinRowsPerProcess = minRowsPerProcess;
		solver.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes, ghosts);
		return solver.NumberOfLevels();
	};
	int perProcess = 1500;
	EXPECT_EQ(levels(1000, perProcess), levels(perProcess * World().Size(), 0));
	EXPECT_EQ(levels(1000, (int)system.Faces->GlobalSize), 1);
}

// The inputs are checked on all the processes: an error on one of them is an error on all of them
TEST(DistributedLibrary, CollectiveErrors)
{
	if (World().Size() == 1)
		GTEST_SKIP() << "several processes";
	int n = 20;
	CartesianHybridSystem global(n, 0);
	PartitionedSystem system(global, n, PartitionShape::Strips, World());
	fhhos4::SparseMatrix A_TT = system.A_TT, A_TF = system.A_TF.LocalColumns();
	fhhos4::SparseMatrix A_FF = DistributedMatrix::RenumberColumns(SparseMatrix(system.Global.A_FF.middleRows(system.Faces->Begin, system.Faces->NOwned)), [&](BigNumber f) { return system.A_TF.ColMap->Local(f); }, system.A_TF.ColMap->NLocal(), 1);
	vector<int64_t> ghosts(system.A_TF.ColMap->Ghosts.begin(), system.A_TF.ColMap->Ghosts.end());
	Vector cellOnes = Vector::Ones(system.Cells->NOwned), faceOnes = Vector::Ones(system.Faces->NOwned);
	// A wrong size of the interpolation of 1 on the last process only
	if (World().Rank() == World().Size() - 1)
		faceOnes = Vector::Ones(system.Faces->NOwned + 1);
	fhhos4::Solver solver = CartesianSolver(World().Comm());
	EXPECT_THROW(solver.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes, ghosts), fhhos4::Error);
	faceOnes = Vector::Ones(system.Faces->NOwned);
	// A ghost beyond 32 bits on the last process, whose low 32 bits are a valid ghost: refused whatever the index type
	vector<int64_t> farGhosts = ghosts;
	if (World().Rank() == World().Size() - 1 && !farGhosts.empty())
		farGhosts[0] += (int64_t)1 << 32;
	EXPECT_THROW(solver.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes, farGhosts), fhhos4::Error);
	// A_FF lower triangular on the last process (as SetupSym() reads it on one process; here diagonal: one coupling
	// added), read in full on several processes
	fhhos4::SparseMatrix lowerA_FF = A_FF;
	if (World().Rank() == World().Size() - 1 && lowerA_FF.rows() > 1)
		lowerA_FF.coeffRef(1, 0) = -0.1;
	lowerA_FF.makeCompressed();
	try
	{
		solver.SetupSym(A_TT, A_TF, lowerA_FF, cellOnes, faceOnes, ghosts);
		ADD_FAILURE() << "a lower triangular A_FF accepted on several processes";
	}
	catch (const fhhos4::Error& e)
	{
		EXPECT_NE(string(e.what()).find("lower triangular"), string::npos) << e.what();
	}
	// Without Comm: each process alone, which refuses the inputs of a distributed system, the ghosts or the columns of
	// ghost faces in A_FF (here those of A_TF: none on rank 0, which owns the faces between the processes), and tells to
	// set Comm
	if (!ghosts.empty())
	{
		fhhos4::Solver noComm = CartesianSolver(MPI_COMM_NULL);
		for (bool withGhosts : { true, false })
		{
			try
			{
				if (withGhosts)
					noComm.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes, ghosts);
				else
					noComm.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes);
				ADD_FAILURE() << "a distributed system accepted without Comm";
			}
			catch (const fhhos4::Error& e)
			{
				EXPECT_NE(string(e.what()).find("set Comm"), string::npos) << e.what();
			}
		}
	}
}
