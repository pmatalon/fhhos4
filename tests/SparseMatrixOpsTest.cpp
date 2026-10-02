// Unit tests of the parallel sparse matrix operations used in the setup of the algebraic multigrids
// (src/Utils/SparseMatrixOps.h): they must give exactly the same matrices as the Eigen operations they replace
// (same structure, same values bit for bit), whatever the number of threads.
#include <gtest/gtest.h>
#include <cstring>
#include <random>
#include <vector>
#include "Utils/SparseMatrixOps.h"

namespace
{
	// Random n x m matrix with about nnzPerRow non-zeros per row, including explicit zeros and empty rows
	SparseMatrix RandomMatrix(BigNumber n, BigNumber m, int nnzPerRow, unsigned seed)
	{
		std::mt19937 gen(seed);
		std::uniform_int_distribution<BigNumber> col(0, m - 1);
		std::uniform_real_distribution<double> value(-1, 1);
		std::vector<Eigen::Triplet<double, SparseMatrixIndex>> triplets;
		for (BigNumber i = 0; i < n; i++)
		{
			if (i % 7 == 3)
				continue;
			for (int k = 0; k < nnzPerRow; k++)
				triplets.push_back({ (SparseMatrixIndex)i, (SparseMatrixIndex)col(gen), k == 0 ? 0.0 : value(gen) });
		}
		SparseMatrix A(n, m);
		A.setFromTriplets(triplets.begin(), triplets.end());
		return A;
	}

	::testing::AssertionResult Identical(const SparseMatrix& X, const SparseMatrix& Y)
	{
		if (!X.isCompressed() || !Y.isCompressed())
			return ::testing::AssertionFailure() << "uncompressed matrix";
		if (X.rows() != Y.rows() || X.cols() != Y.cols() || X.nonZeros() != Y.nonZeros())
			return ::testing::AssertionFailure() << "different sizes or numbers of non-zeros (" << X.nonZeros() << " vs " << Y.nonZeros() << ")";
		if (!std::equal(X.outerIndexPtr(), X.outerIndexPtr() + X.rows() + 1, Y.outerIndexPtr())
			|| !std::equal(X.innerIndexPtr(), X.innerIndexPtr() + X.nonZeros(), Y.innerIndexPtr()))
			return ::testing::AssertionFailure() << "different structures";
		if (std::memcmp(X.valuePtr(), Y.valuePtr(), X.nonZeros() * sizeof(double)) != 0)
			return ::testing::AssertionFailure() << "different values";
		return ::testing::AssertionSuccess();
	}

	const std::vector<int> NThreads = { 1, 3, 8 };
}

TEST(SparseMatrixOpsTest, MultiplyGivesEigenProduct)
{
	SparseMatrix A = RandomMatrix(500, 300, 6, 1);
	SparseMatrix B = RandomMatrix(300, 400, 5, 2);
	SparseMatrix expected = A * B;
	for (int nThreads : NThreads)
	{
		Parallelism::SetNThreads(nThreads);
		EXPECT_TRUE(Identical(SparseMatrixOps::Multiply(A, B), expected)) << nThreads << " threads";
	}
	Parallelism::SetNThreads(0);
}

TEST(SparseMatrixOpsTest, TripleProductGivesEigenProduct)
{
	// Galerkin product as written in the U-AMG setup: P^T * S * P with S given by its lower triangular part
	SparseMatrix S = RandomMatrix(400, 400, 7, 3); // not symmetric: only its lower triangular part must be used
	SparseMatrix P = RandomMatrix(400, 150, 4, 4);
	SparseMatrix expected = P.transpose() * S.selfadjointView<Eigen::Lower>() * P;
	for (int nThreads : NThreads)
	{
		Parallelism::SetNThreads(nThreads);
		SparseMatrix Pt = SparseMatrixOps::Transpose(P);
		SparseMatrix PtSP = SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(Pt, SparseMatrixOps::FullFromLower(S)), P);
		EXPECT_TRUE(Identical(PtSP, expected)) << nThreads << " threads";
	}
	Parallelism::SetNThreads(0);
}

TEST(SparseMatrixOpsTest, FullFromLowerGivesSelfAdjointView)
{
	SparseMatrix S = RandomMatrix(300, 300, 6, 5);
	SparseMatrix expected = S.selfadjointView<Eigen::Lower>();
	expected.makeCompressed();
	for (int nThreads : NThreads)
	{
		Parallelism::SetNThreads(nThreads);
		EXPECT_TRUE(Identical(SparseMatrixOps::FullFromLower(S), expected)) << nThreads << " threads";
	}
	Parallelism::SetNThreads(0);
}

TEST(SparseMatrixOpsTest, SelectRowsGivesCopiedRows)
{
	int blockSize = 3;
	SparseMatrix A = RandomMatrix(300, 200, 5, 6);
	SparseMatrix B = RandomMatrix(300, 200, 5, 7);
	std::vector<bool> fromA(300 / blockSize);
	for (size_t i = 0; i < fromA.size(); i++)
		fromA[i] = (i % 3 != 0);

	// Same rows copied by NonZeroCoefficients (which drops the coefficients <= ZeroThreshold)
	NonZeroCoefficients coeffs;
	for (size_t i = 0; i < fromA.size(); i++)
		coeffs.CopyRows(i * blockSize, blockSize, fromA[i] ? A : B);
	SparseMatrix expected(A.rows(), A.cols());
	coeffs.Fill(expected);

	for (int nThreads : NThreads)
	{
		Parallelism::SetNThreads(nThreads);
		EXPECT_TRUE(Identical(SparseMatrixOps::SelectRows(fromA, A, B, blockSize, NonZeroCoefficients::ZeroThreshold), expected)) << nThreads << " threads";
	}
	Parallelism::SetNThreads(0);
}
