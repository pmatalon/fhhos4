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

// The uncondensed system distributed by rows, through the interface of the library: the ghost faces numbered in reverse
// order of their global numbers (any order). The solution of the system on one process.
TEST(DistributedLibrary, UncondensedSystem)
{
	int n = 60;
	CartesianHybridSystem global(n, 0);
	PartitionedSystem system(global, n, PartitionShape::Blocks, World());
	const IndexMap& faces = *system.Faces;
	const IndexMap& cells = *system.Cells;
	Vector b_T = Vector::LinSpaced(global.A_TT.rows(), 1, 2), b_F = Vector::Zero(global.S.rows());

	// One process: the whole system
	Vector expectedX_T(b_T.rows()), expectedX_F = Vector::Zero(b_F.rows());
	int expectedIterations;
	{
		fhhos4::Solver solver = CartesianSolver(MPI_COMM_SELF);
		solver.SetupSymmetric(system.Global.A_TT, system.Global.A_TF, system.Global.A_FF, Vector::Ones(global.A_TT.rows()), Vector::Ones(global.S.rows()));
		expectedIterations = solver.Solve(b_T, b_F, expectedX_T, expectedX_F).Iterations;
	}

	// The faces of the rows of the process (of A_TF and A_FF): its own, then the others, in reverse order of their global
	// numbers
	vector<BigNumber> ghostFaces;
	SparseMatrix ownA_TF = system.Global.A_TF.middleRows(cells.Begin, cells.NOwned);
	SparseMatrix ownA_FF = system.Global.A_FF.middleRows(faces.Begin, faces.NOwned);
	for (const SparseMatrix* M : { &ownA_TF, &ownA_FF })
	{
		for (BigNumber p = 0; p < (BigNumber)M->nonZeros(); p++)
		{
			BigNumber f = M->innerIndexPtr()[p];
			if (!faces.IsOwned(f))
				ghostFaces.push_back(f);
		}
	}
	sort(ghostFaces.begin(), ghostFaces.end());
	ghostFaces.erase(unique(ghostFaces.begin(), ghostFaces.end()), ghostFaces.end());
	reverse(ghostFaces.begin(), ghostFaces.end());
	auto localFace = [&](BigNumber f) -> BigNumber
	{
		if (faces.IsOwned(f))
			return f - faces.Begin;
		auto it = find(ghostFaces.begin(), ghostFaces.end(), f);
		return it == ghostFaces.end() ? -1 : faces.NOwned + (it - ghostFaces.begin());
	};
	BigNumber nLocalFaces = faces.NOwned + ghostFaces.size();
	auto localColumns = [&](const SparseMatrix& rows)
	{
		return DistributedMatrix::RenumberColumns(rows, localFace, nLocalFaces, 1);
	};
	fhhos4::SparseMatrix A_TT = system.A_TT, A_TF = localColumns(ownA_TF);
	fhhos4::SparseMatrix A_FF = localColumns(ownA_FF);
	vector<int64_t> ghosts(ghostFaces.begin(), ghostFaces.end());

	fhhos4::Solver solver = CartesianSolver(World().Comm());
	solver.SetupSymmetric(A_TT, A_TF, A_FF, Vector::Ones(cells.NOwned), Vector::Ones(faces.NOwned), fhhos4::GhostNumbers(ghosts.size(), ghosts.data()));
	Vector x_T(cells.NOwned), x_F = Vector::Zero(faces.NOwned);
	fhhos4::Result result = solver.Solve(b_T.segment(cells.Begin, cells.NOwned), b_F.segment(faces.Begin, faces.NOwned), x_T, x_F);

	EXPECT_TRUE(result.Converged);
	EXPECT_LE(result.Iterations, expectedIterations + 4);
	if (World().Size() == 1)
		EXPECT_EQ(result.Iterations, expectedIterations);
	Vector allX_T = Gathered(x_T), allX_F = Gathered(x_F);
	EXPECT_LE((allX_F - expectedX_F).norm(), 1e-8 * expectedX_F.norm());
	EXPECT_LE((allX_T - expectedX_T).norm(), 1e-8 * expectedX_T.norm());
	if (World().IsRoot())
		cout << "\t" << World().Size() << " processes: " << result.Iterations << " iterations (" << expectedIterations << " on one)" << endl;
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
	EXPECT_THROW(solver.SetupSymmetric(A_TT, A_TF, A_FF, cellOnes, faceOnes, fhhos4::GhostNumbers(ghosts.size(), ghosts.data())), fhhos4::Error);
	// Without Comm
	fhhos4::Solver noComm = CartesianSolver(MPI_COMM_NULL);
	EXPECT_THROW(noComm.SetupSymmetric(A_TT, A_TF, A_FF, cellOnes, faceOnes, fhhos4::GhostNumbers(ghosts.size(), ghosts.data())), fhhos4::Error);
}
