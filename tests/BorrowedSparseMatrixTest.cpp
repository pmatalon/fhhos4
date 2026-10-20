// Unit tests of BorrowedSparseMatrix (src/Utils/BorrowedSparseMatrix.h): a SparseMatrix that reads the arrays of the
// calling code in place (the inputs of the library fhhos4). It must behave as the matrix of these arrays wherever a
// const SparseMatrix& is read, and leave the arrays untouched when it gives them back or is destroyed (Eigen frees its
// own arrays only: run under valgrind to check the memory).
#include <gtest/gtest.h>
#include <random>
#include <vector>
#include "Utils/BorrowedSparseMatrix.h"
#include "Utils/SparseMatrixOps.h"

namespace
{
	// Random n x m matrix with about nnzPerRow non-zeros per row
	SparseMatrix RandomMatrix(BigNumber n, BigNumber m, int nnzPerRow, unsigned seed)
	{
		std::mt19937 gen(seed);
		std::uniform_int_distribution<BigNumber> col(0, m - 1);
		std::uniform_real_distribution<double> value(-1, 1);
		std::vector<Eigen::Triplet<double, SparseMatrixIndex>> triplets;
		for (BigNumber i = 0; i < n; i++)
			for (int k = 0; k < nnzPerRow; k++)
				triplets.push_back({ (SparseMatrixIndex)i, (SparseMatrixIndex)col(gen), value(gen) });
		SparseMatrix A(n, m);
		A.setFromTriplets(triplets.begin(), triplets.end());
		return A;
	}

	// The arrays of a matrix, as another code holds them
	struct Arrays
	{
		BigNumber Rows, Cols;
		std::vector<SparseMatrixIndex> RowPointers, Columns;
		std::vector<double> Values;

		Arrays(const SparseMatrix& A)
			: Rows(A.rows()), Cols(A.cols()), RowPointers(A.outerIndexPtr(), A.outerIndexPtr() + A.rows() + 1),
			  Columns(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros()), Values(A.valuePtr(), A.valuePtr() + A.nonZeros())
		{}
	};
}

// The borrowed matrix is the matrix of the arrays, read in place, in the products and operations of the solvers
TEST(BorrowedSparseMatrixTest, ReadsTheArraysInPlace)
{
	SparseMatrix A = RandomMatrix(300, 200, 5, 1), B = RandomMatrix(200, 100, 4, 2);
	Arrays a(A);
	BorrowedSparseMatrix borrowed(a.Rows, a.Cols, a.RowPointers.data(), a.Columns.data(), a.Values.data());
	const SparseMatrix& M = borrowed;
	EXPECT_EQ(M.rows(), A.rows());
	EXPECT_EQ(M.cols(), A.cols());
	EXPECT_EQ(M.nonZeros(), A.nonZeros());
	EXPECT_TRUE(M.isCompressed());
	EXPECT_EQ(M.valuePtr(), a.Values.data());
	EXPECT_EQ(M.innerIndexPtr(), a.Columns.data());
	EXPECT_EQ(M.outerIndexPtr(), a.RowPointers.data());

	Eigen::VectorXd x = Eigen::VectorXd::Random(A.cols());
	EXPECT_EQ((M * x - A * x).norm(), 0);
	EXPECT_EQ((SparseMatrixOps::Multiply(M, B) - SparseMatrixOps::Multiply(A, B)).norm(), 0);
	EXPECT_EQ((SparseMatrixOps::Transpose(M) - SparseMatrixOps::Transpose(A)).norm(), 0);
	EXPECT_EQ((SparseMatrix(M) - A).norm(), 0); // a copy: owns its arrays
}

// Borrow() again, GiveBack() and the destruction leave the arrays untouched; the matrix becomes empty again
TEST(BorrowedSparseMatrixTest, GivesTheArraysBack)
{
	SparseMatrix A = RandomMatrix(50, 40, 3, 3), B = RandomMatrix(20, 30, 2, 4);
	Arrays a(A), b(B);
	{
		BorrowedSparseMatrix M;
		EXPECT_EQ(M.rows(), 0);
		M.Borrow(a.Rows, a.Cols, a.RowPointers.data(), a.Columns.data(), a.Values.data());
		EXPECT_EQ(M.valuePtr(), a.Values.data());
		M.Borrow(b.Rows, b.Cols, b.RowPointers.data(), b.Columns.data(), b.Values.data());
		EXPECT_EQ(M.rows(), B.rows());
		EXPECT_EQ((M - B).norm(), 0);
		M.GiveBack();
		EXPECT_EQ(M.rows(), 0);
		EXPECT_EQ(M.cols(), 0);
		EXPECT_EQ(M.nonZeros(), 0);
		M.Borrow(a.Rows, a.Cols, a.RowPointers.data(), a.Columns.data(), a.Values.data());
	} // destroyed while borrowing
	EXPECT_EQ(Arrays(A).Values, a.Values);
	EXPECT_EQ(Arrays(A).Columns, a.Columns);
	EXPECT_EQ(Arrays(A).RowPointers, a.RowPointers);
	EXPECT_EQ(Arrays(B).Values, b.Values);
}

// A matrix without rows (its row pointers: one 0) or without non-zeros
TEST(BorrowedSparseMatrixTest, EmptyMatrices)
{
	SparseMatrixIndex noRows[1] = { 0 };
	BorrowedSparseMatrix empty(0, 5, noRows, nullptr, nullptr);
	EXPECT_EQ(empty.rows(), 0);
	EXPECT_EQ(empty.cols(), 5);
	EXPECT_EQ(empty.nonZeros(), 0);

	std::vector<SparseMatrixIndex> zeros(4, 0);
	BorrowedSparseMatrix noNonZeros(3, 3, zeros.data(), nullptr, nullptr);
	EXPECT_EQ(noNonZeros.nonZeros(), 0);
	EXPECT_EQ((noNonZeros * Eigen::VectorXd::Ones(3)).norm(), 0);
}
