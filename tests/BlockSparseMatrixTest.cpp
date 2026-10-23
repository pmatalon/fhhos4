// Unit tests of the block CSR storage (src/Utils/BlockSparseMatrix.h): built from CSR matrices with blocks of a size
// known at compile time (DispatchBlockSize()) or not, of variable sizes, with coefficients missing in the blocks; its
// product must be that of the CSR matrix summed row by row in the order of the columns (one accumulator, as the CSR
// kernels of the smoothers), bit for bit in double, and to rounding in float.
#include <gtest/gtest.h>
#include <random>
#include <vector>
#include "Utils/BlockSparseMatrix.h"

namespace
{
	// A random sparse matrix on the blocks starting at blockStart: each block row coupled to itself and to a few other
	// blocks (dense blocks), without the coefficients of probability missing (0 by default)
	SparseMatrix BlockMatrix(const std::vector<SparseMatrixIndex>& blockStart, unsigned seed, double missing = 0)
	{
		std::mt19937 gen(seed);
		std::uniform_real_distribution<double> value(-1, 1), uniform(0, 1);
		int nb = blockStart.size() - 1;
		std::uniform_int_distribution<int> otherBlock(0, nb - 1);
		std::vector<Eigen::Triplet<double, SparseMatrixIndex>> triplets;
		for (int I = 0; I < nb; I++)
		{
			std::vector<int> blocks = { I };
			for (int k = 0; k < 4; k++)
				blocks.push_back(otherBlock(gen));
			for (int J : blocks)
				for (SparseMatrixIndex i = blockStart[I]; i < blockStart[I + 1]; i++)
					for (SparseMatrixIndex j = blockStart[J]; j < blockStart[J + 1]; j++)
						if (!(uniform(gen) < missing))
							triplets.push_back({ i, j, value(gen) });
		}
		SparseMatrix A(blockStart.back(), blockStart.back());
		A.setFromTriplets(triplets.begin(), triplets.end()); // the duplicates (blocks drawn twice) summed
		return A;
	}

	std::vector<SparseMatrixIndex> UniformBlocks(int nb, int blockSize)
	{
		std::vector<SparseMatrixIndex> start(nb + 1);
		for (int I = 0; I <= nb; I++)
			start[I] = I * blockSize;
		return start;
	}

	bool Identical(const Vector& x, const Vector& y)
	{
		return x.size() == y.size() && (x.array() == y.array()).all();
	}

	// A x, each row summed in the order of its columns with one accumulator (Eigen 5's product uses two)
	Vector Product(const SparseMatrix& A, const Vector& x)
	{
		Vector y(A.rows());
		for (BigNumber i = 0; i < A.rows(); i++)
		{
			double s = 0;
			for (SparseMatrix::InnerIterator it(A, i); it; ++it)
				s += it.value() * x[it.col()];
			y[i] = s;
		}
		return y;
	}
}

// Blocks of the sizes compiled for (2, 3, 4, 5, 6, 10, 15) and others (1, 7: the sizes at run time): the product of the
// CSR matrix, bit for bit, on 1 and several threads
TEST(BlockSparseMatrixTest, UniformBlocks)
{
	for (int blockSize : { 1, 2, 3, 4, 5, 6, 7, 10, 15 })
	{
		SparseMatrix A = BlockMatrix(UniformBlocks(1500, blockSize), blockSize);
		Vector x = Vector::Random(A.cols());
		Vector expected = Product(A, x);
		BlockSparseMatrix<double> B = BlockSparseMatrix<double>::FromCSR(A, blockSize);
		EXPECT_EQ(B.BlockSize, blockSize);
		EXPECT_EQ(B.Rows(), A.rows());
		EXPECT_EQ((BigNumber)B.NonZeros(), (BigNumber)A.nonZeros()) << "full blocks: no zero added";
		for (int threads : { 1, 4 })
		{
			Parallelism::SetNThreads(threads);
			EXPECT_TRUE(Identical(B * x, expected)) << "block size " << blockSize << ", " << threads << " threads";
		}
		Parallelism::SetNThreads(0);
	}
}

// Blocks of different sizes (variable degrees): the sizes at run time, the same product
TEST(BlockSparseMatrixTest, VariableBlocks)
{
	std::mt19937 gen(3);
	std::uniform_int_distribution<int> size(1, 6);
	std::vector<SparseMatrixIndex> start = { 0 };
	for (int I = 0; I < 2000; I++)
		start.push_back(start.back() + size(gen));
	SparseMatrix A = BlockMatrix(start, 4);
	Vector x = Vector::Random(A.cols());
	BlockSparseMatrix<double> B = BlockSparseMatrix<double>::FromCSR(A, start);
	EXPECT_EQ(B.BlockSize, 0);
	EXPECT_EQ((BigNumber)B.NonZeros(), (BigNumber)A.nonZeros());
	EXPECT_TRUE(Identical(B * x, Product(A, x)));
	// The diagonal blocks
	for (BigNumber I : { 0, 7, 1999 })
		EXPECT_EQ(B.DiagonalBlock(I), DenseMatrix(A.block(start[I], start[I], start[I + 1] - start[I], start[I + 1] - start[I])));
}

// Coefficients missing in the blocks: stored as zeros, the same product
TEST(BlockSparseMatrixTest, MissingCoefficients)
{
	SparseMatrix A = BlockMatrix(UniformBlocks(1000, 3), 5, 0.3);
	BlockSparseMatrix<double> B = BlockSparseMatrix<double>::FromCSR(A, 3);
	EXPECT_GT((BigNumber)B.NonZeros(), (BigNumber)A.nonZeros());
	Vector x = Vector::Random(A.cols());
	EXPECT_TRUE(Identical(B * x, Product(A, x)));
}

// In single precision: the values rounded, the sums in double
TEST(BlockSparseMatrixTest, SinglePrecision)
{
	SparseMatrix A = BlockMatrix(UniformBlocks(1500, 6), 6);
	BlockSparseMatrix<float> B = BlockSparseMatrix<float>::FromCSR(A, 6);
	FloatVector x = FloatVector::Random(A.cols());
	Vector expected = A * x.cast<double>();
	FloatVector y = B * x;
	EXPECT_LT((y.cast<double>() - expected).norm() / expected.norm(), 1e-6);
}
