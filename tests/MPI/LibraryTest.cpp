#include <gtest/gtest.h>
#include <fhhos4>
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

// The uncondensed system distributed by rows, through the interface of the library: the ghost faces and cells numbered
// in reverse order of their global numbers (any order). The solution of the system on one process. Non-symmetric if
// asymmetry != 0. singular: the pure Neumann system (SingularMatrix; the right-hand side incompatible: projected),
// the solutions compared up to the kernel (the constants, of each component). components: the system on each of them
// (the blocks ⊗ I_n, VectorComponents = n; the ghosts are the same faces and cells).
void ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation prolongation, double asymmetry = 0, bool singular = false, int components = 1)
{
	int n = 60;
	CartesianHybridSystem global(n, asymmetry, singular);
	PartitionedSystem system(global, n, PartitionShape::Blocks, World());
	const IndexMap& faces = *system.Faces;
	const IndexMap& cells = *system.Cells;
	int c = components;
	Vector b_T = OnEachComponent(Vector::LinSpaced(global.A_TT.rows(), 1, 2), c), b_F = Vector::Zero(c * global.S.rows());

	// One process: the whole system
	Vector expectedX_T(b_T.rows()), expectedX_F = Vector::Zero(b_F.rows());
	int expectedIterations;
	{
		fhhos4::Solver solver = CartesianSolver(MPI_COMM_SELF);
		solver.Prolongation = prolongation;
		solver.SingularMatrix = singular;
		solver.VectorComponents = c;
		Vector cellOnes = Vector::Ones(global.A_TT.rows()), faceOnes = Vector::Ones(global.S.rows());
		fhhos4::SparseMatrix A_TT = OnEachComponent(system.Global.A_TT, c), A_TF = OnEachComponent(system.Global.A_TF, c);
		fhhos4::SparseMatrix A_FT = OnEachComponent(system.Global.A_FT, c), A_FF = OnEachComponent(system.Global.A_FF, c);
		if (asymmetry == 0)
			solver.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes);
		else
			solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellOnes, faceOnes);
		expectedIterations = solver.Solve(b_T, b_F, expectedX_T, expectedX_F).Iterations;
	}

	// The faces of the rows of the process (of A_TF and A_FF): its own, then the others; the cells of A_FT: the same
	SparseMatrix ownA_TF = system.Global.A_TF.middleRows(cells.Begin, cells.NOwned);
	SparseMatrix ownA_FF = system.Global.A_FF.middleRows(faces.Begin, faces.NOwned);
	SparseMatrix ownA_FT = system.Global.A_FT.middleRows(faces.Begin, faces.NOwned);
	vector<BigNumber> ghostFaces = ReversedGhosts({ &ownA_TF, &ownA_FF }, faces), ghostCells = ReversedGhosts({ &ownA_FT }, cells);
	fhhos4::SparseMatrix A_TT = OnEachComponent(system.A_TT, c), A_TF = OnEachComponent(LocalColumns(ownA_TF, faces, ghostFaces), c);
	fhhos4::SparseMatrix A_FF = OnEachComponent(LocalColumns(ownA_FF, faces, ghostFaces), c), A_FT = OnEachComponent(LocalColumns(ownA_FT, cells, ghostCells), c);
	vector<int64_t> faceGhosts(ghostFaces.begin(), ghostFaces.end()), cellGhosts(ghostCells.begin(), ghostCells.end());

	fhhos4::Solver solver = CartesianSolver(World().Comm());
	solver.Prolongation = prolongation;
	solver.SingularMatrix = singular;
	solver.VectorComponents = c;
	Vector cellOnes = Vector::Ones(cells.NOwned), faceOnes = Vector::Ones(faces.NOwned);
	if (asymmetry == 0)
		solver.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes, fhhos4::GhostNumbers(faceGhosts.size(), faceGhosts.data()));
	else
		solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellOnes, faceOnes, fhhos4::GhostNumbers(faceGhosts.size(), faceGhosts.data()), fhhos4::GhostNumbers(cellGhosts.size(), cellGhosts.data()));
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
	if (World().IsRoot())
		cout << "\t" << World().Size() << " processes: " << result.Iterations << " iterations (" << expectedIterations << " on one)" << endl;
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
	EXPECT_THROW(solver.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes, fhhos4::GhostNumbers(ghosts.size(), ghosts.data())), fhhos4::Error);
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
					noComm.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes, fhhos4::GhostNumbers(ghosts.size(), ghosts.data()));
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
