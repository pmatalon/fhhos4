// Unit tests of the reductions of the solvers by blocks of rows (src/Utils/Communicator.h, namespace Parallel): the same
// results bit for bit whatever the number of threads; on one block (at most 4096 rows), Eigen's whole-vector ones; a
// pass that updates then reduces each block (the fused passes of the Krylov methods) gives the separate operations'.
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <vector>
#include "Utils/Communicator.h"

namespace
{
	Vector RandomVector(BigNumber n, unsigned seed)
	{
		std::mt19937 gen(seed);
		std::uniform_real_distribution<double> value(-1, 1);
		Vector x(n);
		for (BigNumber i = 0; i < n; i++)
			x[i] = value(gen);
		return x;
	}

	const std::vector<int> NThreads = { 1, 3, 8 };
	// One block, several, a last partial block, and blocks larger than 4096 rows (above 4M rows)
	const std::vector<BigNumber> Sizes = { 0, 1, 4096, 4097, 100000, 4096 * 1024 + 5 };
}

TEST(ParallelReductionsTest, SameResultsWhateverTheNumberOfThreads)
{
	for (BigNumber n : Sizes)
	{
		Vector x = RandomVector(n, 1);
		Vector y = RandomVector(n, 2);
		FloatVector xf = x.cast<float>();
		FloatVector yf = y.cast<float>();
		Parallelism::SetNThreads(1);
		double dot = Parallel::Dot(x, y);
		double squaredNorm = Parallel::SquaredNorm(x);
		double dotFloat = Parallel::Dot(xf, yf);
		for (int nThreads : NThreads)
		{
			Parallelism::SetNThreads(nThreads);
			EXPECT_EQ(Parallel::Dot(x, y), dot) << n << " rows, " << nThreads << " threads";
			EXPECT_EQ(Parallel::SquaredNorm(x), squaredNorm) << n << " rows, " << nThreads << " threads";
			EXPECT_EQ(Parallel::Norm(x), std::sqrt(squaredNorm)) << n << " rows, " << nThreads << " threads";
			EXPECT_EQ(Parallel::Dot(xf, yf), dotFloat) << n << " rows, " << nThreads << " threads";
		}
		if (n <= 4096)
		{
			EXPECT_EQ(dot, x.dot(y)) << n << " rows";
			EXPECT_EQ(squaredNorm, x.squaredNorm()) << n << " rows";
		}
		EXPECT_NEAR(dot, x.dot(y), 1e-12 * std::max<double>(1, n)) << n << " rows";
	}
	Parallelism::SetNThreads(0);
}

TEST(ParallelReductionsTest, RowBlocks)
{
	for (BigNumber n : Sizes)
	{
		Parallel::RowBlocks blocks(n);
		EXPECT_LE(blocks.Count, Parallel::RowBlocks::MaxCount) << n << " rows";
		EXPECT_GE(blocks.Count * blocks.Size, n) << n << " rows";
		EXPECT_LT((blocks.Count - 1) * blocks.Size, std::max<BigNumber>(n, 1)) << n << " rows";
		EXPECT_EQ(blocks.Size % 4096, 0) << n << " rows";
	}
	EXPECT_EQ(Parallel::RowBlocks(4096 * 1024).Size, 4096);
	EXPECT_EQ(Parallel::RowBlocks(4096 * 1024 + 1).Size, 8192);
}

TEST(ParallelReductionsTest, FusedPassGivesSeparateOperations)
{
	BigNumber n = 100000;
	Vector x = RandomVector(n, 3);
	Vector y = RandomVector(n, 4);
	double a = 0.37;
	// Separate: y += a x, then y.x and |y|^2
	Vector ySeparate = y;
	ySeparate += a * x;
	double dot = Parallel::Dot(ySeparate, x);
	double squaredNorm = Parallel::SquaredNorm(ySeparate);
	for (int nThreads : NThreads)
	{
		Parallelism::SetNThreads(nThreads);
		Vector yFused = y;
		std::array<double, 2> sums = Parallel::SumsOverRowBlocks<2>(n, [&](BigNumber i, BigNumber m)
		{
			yFused.segment(i, m) += a * x.segment(i, m);
			return std::array<double, 2>{ Parallel::DotOfRows(yFused, x, i, m), Parallel::SquaredNormOfRows(yFused, i, m) };
		});
		EXPECT_TRUE(yFused == ySeparate) << nThreads << " threads";
		EXPECT_EQ(sums[0], dot) << nThreads << " threads";
		EXPECT_EQ(sums[1], squaredNorm) << nThreads << " threads";
	}
	Parallelism::SetNThreads(0);
}
