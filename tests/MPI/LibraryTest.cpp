#include <gtest/gtest.h>
#include <fhhos4>
#include "support/CartesianHybridSystem.h"
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
// the solutions compared up to the kernel (the constants).
void ExpectSolvesTheUncondensedSystem(fhhos4::Prolongation prolongation, double asymmetry = 0, bool singular = false)
{
	int n = 60;
	CartesianHybridSystem global(n, asymmetry, singular);
	PartitionedSystem system(global, n, PartitionShape::Blocks, World());
	const IndexMap& faces = *system.Faces;
	const IndexMap& cells = *system.Cells;
	Vector b_T = Vector::LinSpaced(global.A_TT.rows(), 1, 2), b_F = Vector::Zero(global.S.rows());

	// One process: the whole system
	Vector expectedX_T(b_T.rows()), expectedX_F = Vector::Zero(b_F.rows());
	int expectedIterations;
	{
		fhhos4::Solver solver = CartesianSolver(MPI_COMM_SELF);
		solver.Prolongation = prolongation;
		solver.SingularMatrix = singular;
		Vector cellOnes = Vector::Ones(global.A_TT.rows()), faceOnes = Vector::Ones(global.S.rows());
		if (asymmetry == 0)
			solver.SetupSym(system.Global.A_TT, system.Global.A_TF, system.Global.A_FF, cellOnes, faceOnes);
		else
			solver.SetupNonSym(system.Global.A_TT, system.Global.A_TF, system.Global.A_FT, system.Global.A_FF, cellOnes, faceOnes);
		expectedIterations = solver.Solve(b_T, b_F, expectedX_T, expectedX_F).Iterations;
	}

	// The faces of the rows of the process (of A_TF and A_FF): its own, then the others; the cells of A_FT: the same
	SparseMatrix ownA_TF = system.Global.A_TF.middleRows(cells.Begin, cells.NOwned);
	SparseMatrix ownA_FF = system.Global.A_FF.middleRows(faces.Begin, faces.NOwned);
	SparseMatrix ownA_FT = system.Global.A_FT.middleRows(faces.Begin, faces.NOwned);
	vector<BigNumber> ghostFaces = ReversedGhosts({ &ownA_TF, &ownA_FF }, faces), ghostCells = ReversedGhosts({ &ownA_FT }, cells);
	fhhos4::SparseMatrix A_TT = system.A_TT, A_TF = LocalColumns(ownA_TF, faces, ghostFaces);
	fhhos4::SparseMatrix A_FF = LocalColumns(ownA_FF, faces, ghostFaces), A_FT = LocalColumns(ownA_FT, cells, ghostCells);
	vector<int64_t> faceGhosts(ghostFaces.begin(), ghostFaces.end()), cellGhosts(ghostCells.begin(), ghostCells.end());

	fhhos4::Solver solver = CartesianSolver(World().Comm());
	solver.Prolongation = prolongation;
	solver.SingularMatrix = singular;
	Vector cellOnes = Vector::Ones(cells.NOwned), faceOnes = Vector::Ones(faces.NOwned);
	if (asymmetry == 0)
		solver.SetupSym(A_TT, A_TF, A_FF, cellOnes, faceOnes, fhhos4::GhostNumbers(faceGhosts.size(), faceGhosts.data()));
	else
		solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellOnes, faceOnes, fhhos4::GhostNumbers(faceGhosts.size(), faceGhosts.data()), fhhos4::GhostNumbers(cellGhosts.size(), cellGhosts.data()));
	Vector x_T(cells.NOwned), x_F = Vector::Zero(faces.NOwned);
	fhhos4::Result result = solver.Solve(b_T.segment(cells.Begin, cells.NOwned), b_F.segment(faces.Begin, faces.NOwned), x_T, x_F);

	EXPECT_TRUE(result.Converged);
	// singular: with HarmonicExtension, 14 iterations on one process (and on 3), 18-19 on 2 and 4, as the regular system
	EXPECT_LE(result.Iterations, expectedIterations + (singular ? 6 : 4));
	if (World().Size() == 1)
		EXPECT_EQ(result.Iterations, expectedIterations);
	Vector allX_T = Gathered(x_T), allX_F = Gathered(x_F);
	if (singular)
	{
		// The difference without its component along the kernel (the interpolation of 1: ones)
		double shift = ((allX_T - expectedX_T).sum() + (allX_F - expectedX_F).sum()) / (allX_T.rows() + allX_F.rows());
		allX_T.array() -= shift;
		allX_F.array() -= shift;
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
