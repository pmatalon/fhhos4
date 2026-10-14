// Unit tests of the static condensation by cells (src/Solver/StaticCondensation.h), used by the library fhhos4
// to condense the uncondensed systems it is given: the condensed matrix must be bit for bit that of the global products
// FullFromLower(A_FF) - Multiply(A_FT, Multiply(inv(A_TT), A_TF)) (inverses of the blocks in full), whatever the number
// of threads, in the pattern of A_FF (one pass) or not (global products); the vectors must solve the uncondensed system.
#include <gtest/gtest.h>
#include <cstring>
#include <random>
#include <vector>
#include "Solver/StaticCondensation.h"

namespace
{
	// Hybrid system of nCells cells in a row, each with 3 faces: the face shared with the previous cell, the one shared
	// with the next cell, and one of its own. cellBlockSize DoFs per cell, faceBlockSize per face. The local matrix of a
	// cell is dense, positive definite (its symmetric part), with a skew-symmetric part if !symmetric. A_FF has the
	// couplings of the faces of each cell (as HHO), or only its diagonal blocks (faceCouplings false).
	struct BlockSystem
	{
		int CellBlockSize, FaceBlockSize;
		SparseMatrix A_TT, A_TF, A_FT, A_FF;

		BlockSystem(int nCells, int cellBlockSize, int faceBlockSize, bool symmetric, bool faceCouplings, unsigned seed)
			: CellBlockSize(cellBlockSize), FaceBlockSize(faceBlockSize)
		{
			int nT = cellBlockSize, nF = faceBlockSize, nFaces = 2 * nCells + 1;
			std::mt19937 gen(seed);
			std::uniform_real_distribution<double> value(-1, 1);
			std::vector<Eigen::Triplet<double, SparseMatrixIndex>> TT, TF, FT, FF;
			for (int T = 0; T < nCells; T++)
			{
				int faces[3] = { T, T + 1, nCells + 1 + T };
				int n = nT + 3 * nF;
				DenseMatrix B(n, n), S(n, n);
				for (int i = 0; i < n; i++)
					for (int j = 0; j < n; j++)
					{
						B(i, j) = value(gen);
						S(i, j) = value(gen);
					}
				DenseMatrix local = B * B.transpose() + n * DenseMatrix::Identity(n, n);
				if (!symmetric)
					local += S - S.transpose();
				// Local DoF of the cell or face DoF k of the local face f
				auto global = [&](int f, int k) { return faces[f] * nF + k; };
				for (int i = 0; i < nT; i++)
					for (int j = 0; j < nT; j++)
						TT.push_back({ T * nT + i, T * nT + j, local(i, j) });
				for (int f = 0; f < 3; f++)
					for (int k = 0; k < nF; k++)
						for (int i = 0; i < nT; i++)
						{
							TF.push_back({ T * nT + i, global(f, k), local(i, nT + f * nF + k) });
							FT.push_back({ global(f, k), T * nT + i, local(nT + f * nF + k, i) });
						}
				for (int f = 0; f < 3; f++)
					for (int g = 0; g < 3; g++)
						if (faceCouplings || f == g)
							for (int k = 0; k < nF; k++)
								for (int l = 0; l < nF; l++)
									FF.push_back({ global(f, k), global(g, l), local(nT + f * nF + k, nT + g * nF + l) });
			}
			A_TT = SparseMatrix(nCells * nT, nCells * nT);
			A_TT.setFromTriplets(TT.begin(), TT.end());
			A_TF = SparseMatrix(nCells * nT, nFaces * nF);
			A_TF.setFromTriplets(TF.begin(), TF.end());
			A_FT = SparseMatrix(nFaces * nF, nCells * nT);
			A_FT.setFromTriplets(FT.begin(), FT.end());
			A_FF = SparseMatrix(nFaces * nF, nFaces * nF);
			A_FF.setFromTriplets(FF.begin(), FF.end());
			if (symmetric) // given as lower triangular parts, read so
			{
				A_TT = SparseMatrix(A_TT.triangularView<Eigen::Lower>());
				A_FF = SparseMatrix(A_FF.triangularView<Eigen::Lower>());
			}
		}

		// A_FT A_TT^{-1} A_TF with the global products, the inverses of the blocks in full
		SparseMatrix CellElimination(bool symmetric) const
		{
			SparseMatrix fullA_TT = symmetric ? SparseMatrixOps::FullFromLower(A_TT) : A_TT;
			std::vector<Eigen::Triplet<double, SparseMatrixIndex>> t;
			int nT = CellBlockSize;
			for (BigNumber T = 0; T < (BigNumber)A_TT.rows() / nT; T++)
			{
				DenseMatrix block = fullA_TT.block(T * nT, T * nT, nT, nT);
				DenseMatrix inv = block.inverse();
				for (int i = 0; i < nT; i++)
					for (int j = 0; j < nT; j++)
						t.push_back({ (SparseMatrixIndex)(T * nT + i), (SparseMatrixIndex)(T * nT + j), inv(i, j) });
			}
			SparseMatrix invA_TT(A_TT.rows(), A_TT.cols());
			invA_TT.setFromTriplets(t.begin(), t.end());
			SparseMatrix A_FTsym = symmetric ? SparseMatrixOps::Transpose(A_TF) : A_FT;
			return SparseMatrixOps::Multiply(A_FTsym, SparseMatrixOps::Multiply(invA_TT, A_TF));
		}

		SparseMatrix CondensedMatrix(bool symmetric) const
		{
			SparseMatrix A = (symmetric ? SparseMatrixOps::FullFromLower(A_FF) : A_FF) - CellElimination(symmetric);
			A.makeCompressed();
			return A;
		}
	};

	::testing::AssertionResult Identical(const SparseMatrix& X, const SparseMatrix& Y)
	{
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

// A_FF with the couplings of the faces of each cell: the rows of A in its pattern, in one pass
TEST(StaticCondensationTest, CondensedMatrixInThePatternOfTheFaceBlock)
{
	for (bool symmetric : { true, false })
	{
		BlockSystem system(300, 4, 3, symmetric, true, 1);
		SparseMatrix expected = system.CondensedMatrix(symmetric);
		for (int nThreads : NThreads)
		{
			Parallelism::SetNThreads(nThreads);
			StaticCondensation condensation(system.A_TT, system.A_TF, symmetric ? nullptr : &system.A_FT, system.CellBlockSize);
			EXPECT_TRUE(Identical(condensation.CondensedMatrix(system.A_FF, symmetric), expected)) << (symmetric ? "symmetric, " : "non-symmetric, ") << nThreads << " threads";
		}
	}
	Parallelism::SetNThreads(0);
}

// A_FF without the couplings of the faces of each cell (block diagonal): the global products
TEST(StaticCondensationTest, CondensedMatrixOutsideThePatternOfTheFaceBlock)
{
	for (bool symmetric : { true, false })
	{
		BlockSystem system(300, 3, 2, symmetric, false, 2);
		SparseMatrix expected = system.CondensedMatrix(symmetric);
		for (int nThreads : NThreads)
		{
			Parallelism::SetNThreads(nThreads);
			StaticCondensation condensation(system.A_TT, system.A_TF, symmetric ? nullptr : &system.A_FT, system.CellBlockSize);
			EXPECT_TRUE(Identical(condensation.CondensedMatrix(system.A_FF, symmetric), expected)) << (symmetric ? "symmetric, " : "non-symmetric, ") << nThreads << " threads";
		}
	}
	Parallelism::SetNThreads(0);
}

// A_FF = A + A_FT A_TT^{-1} A_TF, from the condensed matrix A (read in full)
TEST(StaticCondensationTest, FaceBlock)
{
	for (bool symmetric : { true, false })
	{
		BlockSystem system(300, 4, 3, symmetric, true, 3);
		StaticCondensation condensation(system.A_TT, system.A_TF, symmetric ? nullptr : &system.A_FT, system.CellBlockSize);
		SparseMatrix A = condensation.CondensedMatrix(system.A_FF, symmetric);
		SparseMatrix expected = A + system.CellElimination(symmetric);
		expected.makeCompressed();
		EXPECT_TRUE(Identical(condensation.FaceBlock(A), expected)) << (symmetric ? "symmetric" : "non-symmetric");
	}
}

// The right-hand side of the condensed system and the cell unknowns solve the uncondensed system
TEST(StaticCondensationTest, SolvesTheUncondensedSystem)
{
	for (bool symmetric : { true, false })
	{
		BlockSystem system(300, 4, 3, symmetric, true, 4);
		StaticCondensation condensation(system.A_TT, system.A_TF, symmetric ? nullptr : &system.A_FT, system.CellBlockSize);
		SparseMatrix A = condensation.CondensedMatrix(system.A_FF, symmetric);
		SparseMatrix A_TT = symmetric ? SparseMatrixOps::FullFromLower(system.A_TT) : system.A_TT;
		SparseMatrix A_FF = symmetric ? SparseMatrixOps::FullFromLower(system.A_FF) : system.A_FF;
		SparseMatrix A_FT = symmetric ? SparseMatrixOps::Transpose(system.A_TF) : system.A_FT;

		// The solution of the uncondensed system for a given x_F and b_T
		Vector x_F = Vector::Random(A.rows()), b_T = Vector::Random(A_TT.rows());
		Vector invA_TT_b_T = condensation.SolveCells(b_T);
		Vector x_T = condensation.CellUnknowns(invA_TT_b_T, x_F);
		EXPECT_LT((A_TT * x_T + system.A_TF * x_F - b_T).norm() / b_T.norm(), 1e-13);
		Vector b_F = A_FT * x_T + A_FF * x_F;
		Vector b = condensation.CondensedRHS(b_F, invA_TT_b_T);
		EXPECT_LT((A * x_F - b).norm() / b.norm(), 1e-13);
	}
}
