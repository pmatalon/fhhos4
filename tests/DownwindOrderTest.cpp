// Unit tests of the downwind order of the smoothers of U-AMG for a non-symmetric matrix (src/Solver/Multigrid/
// DownwindOrder.h), and of the Gauss-Seidel sweeps in a given order (GaussSeidel, BlockGaussSeidel: SweepOrder).
#include <gtest/gtest.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include "Solver/Multigrid/DownwindOrder.h"
#include "Solver/FixedPoint/BlockGaussSeidel.h"
#include "Solver/FixedPoint/GaussSeidel.h"

namespace
{
	// Upwind finite differences of -eps u'' + beta u' on n points (Dirichlet), in blocks of blockSize DoFs whose first
	// DoF carries the scheme (the others: identity): row i couples to i-1 with -(eps + beta) and to i+1 with -eps for
	// beta > 0 (flow towards increasing i), symmetrically for beta < 0. periodic: the first and last points are neighbours
	// (circulating flow). eps: one coefficient per interval (heterogeneous diffusion).
	SparseMatrix Upwind1D(int n, const std::vector<double>& eps, double beta, bool periodic = false, int blockSize = 1)
	{
		std::vector<Eigen::Triplet<double, SparseMatrixIndex>> triplets;
		auto add = [&](int i, int j, double v) { triplets.push_back({ (SparseMatrixIndex)(i * blockSize), (SparseMatrixIndex)(j * blockSize), v }); };
		for (int i = 0; i < n; i++)
		{
			int left = periodic ? (i + n - 1) % n : i - 1, right = periodic ? (i + 1) % n : i + 1;
			double epsLeft = eps[i % eps.size()], epsRight = eps[(i + 1) % eps.size()];
			double diag = epsLeft + epsRight + std::abs(beta);
			add(i, i, diag);
			if (left >= 0)
				add(i, left, -epsLeft - std::max(beta, 0.0));
			if (right < n)
				add(i, right, -epsRight - std::max(-beta, 0.0));
			for (int k = 1; k < blockSize; k++)
				triplets.push_back({ (SparseMatrixIndex)(i * blockSize + k), (SparseMatrixIndex)(i * blockSize + k), 1.0 });
		}
		SparseMatrix A(n * blockSize, n * blockSize);
		A.setFromTriplets(triplets.begin(), triplets.end());
		return A;
	}

	std::vector<BigNumber> Order(const SparseMatrix& A, int blockSize = 1, const Vector& interpOfOne = Vector())
	{
		BigNumber nOneSided, nCycleBreaks;
		return DownwindOrder(A, blockSize, interpOfOne, 0.25, nOneSided, nCycleBreaks);
	}

	std::vector<BigNumber> Range(BigNumber n, bool reversed = false)
	{
		std::vector<BigNumber> r(n);
		std::iota(r.begin(), r.end(), 0);
		if (reversed)
			std::reverse(r.begin(), r.end());
		return r;
	}
}

// Convection-dominated: the order of the flow, in either direction
TEST(DownwindOrderTest, FollowsTheFlow)
{
	EXPECT_EQ(Order(Upwind1D(10, { 1 }, 10)), Range(10));
	EXPECT_EQ(Order(Upwind1D(10, { 1 }, -10)), Range(10, true));
	// Blocks: the coupling of their first DoFs
	EXPECT_EQ(Order(Upwind1D(10, { 1 }, -10, false, 3), 3), Range(10, true));
}

// No one-sided coupling: a symmetric matrix, even when its rows have different maxima (heterogeneous diffusion, whose
// strong couplings in Notay's sense are not symmetric), or a diffusion-dominated one (cell Peclet number < 1/3)
TEST(DownwindOrderTest, NoneWithoutConvection)
{
	EXPECT_TRUE(Order(Upwind1D(10, { 1, 100, 1, 1, 1000 }, 0)).empty());
	EXPECT_TRUE(Order(Upwind1D(10, { 1 }, 0.1)).empty());
}

// The scaling of the constant modes (interpolation of 1) does not change the order
TEST(DownwindOrderTest, InvariantToTheScalingOfTheBases)
{
	SparseMatrix A = Upwind1D(10, { 1 }, 10);
	Vector c(10);
	for (int i = 0; i < 10; i++)
		c[i] = std::pow(10, (i % 3) - 1);
	SparseMatrix scaled = c.cwiseInverse().asDiagonal() * A * c.cwiseInverse().asDiagonal(); // coordinates z = y/c
	EXPECT_EQ(Order(scaled, 1, c), Range(10));
}

// Circulating flow: one cycle, broken once
TEST(DownwindOrderTest, BreaksTheCycles)
{
	BigNumber nOneSided, nCycleBreaks;
	std::vector<BigNumber> order = DownwindOrder(Upwind1D(10, { 1 }, 10, true), 1, Vector(), 0.25, nOneSided, nCycleBreaks);
	EXPECT_EQ(order, Range(10));
	EXPECT_EQ(nOneSided, 10);
	EXPECT_EQ(nCycleBreaks, 1);
}

// Pure convection (triangular in the downwind order): one Gauss-Seidel sweep in that order solves the system, forward or
// backward (all the sweeps follow the order), sequential or hybrid, scalar or by blocks
TEST(DownwindOrderTest, OneSweepSolvesPureConvection)
{
	for (int blockSize : { 1, 3 })
	{
		SparseMatrix A = Upwind1D(50, { 0 }, -1, false, blockSize); // flow towards decreasing i: upper triangular
		Vector exact = Vector::LinSpaced(A.rows(), 1, 2);
		Vector b = A * exact;
		std::vector<BigNumber> order = Order(A, blockSize);
		ASSERT_EQ(order, Range(50, true));
		for (Direction direction : { Direction::Forward, Direction::Backward })
		{
			for (bool hybrid : { false, true })
			{
				IterativeSolver* gs;
				if (blockSize == 1)
				{
					GaussSeidel* scalar = new GaussSeidel(direction, hybrid);
					scalar->SweepOrder = &order;
					gs = scalar;
				}
				else
				{
					BlockGaussSeidel* block = new BlockGaussSeidel(blockSize, direction, hybrid);
					block->SweepOrder = &order;
					gs = block;
				}
				gs->StoppingCrit = StoppingCriteria::MaxIterations;
				gs->MaxIterations = 1;
				gs->PrintIterationResults = false;
				gs->Setup(A);
				Vector x = Vector::Zero(A.rows());
				gs->Solve(b, x, true);
				// Hybrid: each thread sweeps its chunk of the order: the first block row of the order is exact
				double error = hybrid ? (x - exact).segment(order[0] * blockSize, blockSize).norm() : (x - exact).norm();
				EXPECT_LT(error, 1e-12 * exact.norm()) << "block size " << blockSize << ", hybrid " << hybrid;
				EXPECT_FALSE(gs->CanOptimizeResidualComputation());
				delete gs;
			}
		}
	}
}
