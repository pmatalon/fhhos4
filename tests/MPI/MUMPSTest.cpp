#include <gtest/gtest.h>
#include <Eigen/SparseLU>
#include "Solver/Direct/MUMPS.h"

namespace
{
	const Communicator& World()
	{
		return Communicator::Current();
	}

	// A tridiagonal matrix of size n, symmetric positive definite or not, whose values differ
	SparseMatrix Tridiagonal(BigNumber n, bool spd)
	{
		vector<Eigen::Triplet<double, SparseMatrixIndex>> t;
		for (BigNumber i = 0; i < n; i++)
		{
			t.emplace_back(i, i, 2.5 + 0.01 * i);
			if (i > 0)
				t.emplace_back(i, i - 1, -1.0 - 0.001 * (spd ? i - 1 : i));
			if (i < n - 1)
				t.emplace_back(i, i + 1, spd ? -1.0 - 0.001 * i : -1.0 + 0.002 * i);
		}
		SparseMatrix A(n, n);
		A.setFromTriplets(t.begin(), t.end());
		return A;
	}

	Vector ExactSolution(const SparseMatrix& A, const Vector& b)
	{
		ColMajorSparseMatrix colMajor = A;
		Eigen::SparseLU<ColMajorSparseMatrix> lu(colMajor);
		return lu.solve(b);
	}
}

// MUMPS on the rows of each process (none on the first one, MUMPS's host, on several processes), with blocks of 1 and 2
// rows: each process gets its rows of the solution, at each solve
TEST(MUMPS, DistributedSolve)
{
	for (bool spd : { true, false })
	{
		for (int bs : { 1, 2 })
		{
			BigNumber nOwned = World().IsDistributed() && World().IsRoot() ? 0 : 20 + 7 * World().Rank();
			IndexMap owners(World(), nOwned);
			auto rowMap = make_shared<const IndexMap>(owners);
			SparseMatrix A = Tridiagonal(owners.GlobalSize * bs, spd);
			BigNumber first = owners.Begin * bs, n = nOwned * bs;
			DistributedMatrix dA = DistributedMatrix::FromGlobalColumns(rowMap, owners, A.middleRows(first, n), bs, bs);
			MUMPSSolver solver(spd);
			solver.Setup(dA);
			Vector b = Vector::LinSpaced(A.rows(), 1.0, 2.0);
			Vector exact = ExactSolution(A, b);
			for (double factor : { 1.0, -3.0 }) // the second solve reuses the exchange of the solution
			{
				Vector x = solver.Solve(factor * b.segment(first, n));
				EXPECT_LE((x - factor * exact.segment(first, n)).norm(), 1e-12 * exact.norm()) << "spd " << spd << ", bs " << bs;
			}
		}
	}
}

// The whole matrix on each process: MUMPS on this process alone
TEST(MUMPS, SequentialSolve)
{
	for (bool spd : { true, false })
	{
		SparseMatrix A = Tridiagonal(50, spd);
		MUMPSSolver solver(spd);
		solver.Setup(A);
		Vector b = Vector::LinSpaced(A.rows(), 1.0, 2.0);
		Vector exact = ExactSolution(A, b);
		EXPECT_LE((solver.Solve(b) - exact).norm(), 1e-12 * exact.norm()) << "spd " << spd;
	}
}
