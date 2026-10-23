#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include "Types.h"
#include "Parallelism.h"
using namespace std;

// Operations on row-major sparse matrices, computed in parallel by rows (Eigen's sparse products are sequential).
// Each row is computed independently of the others, so the results do not depend on the number of threads.
namespace SparseMatrixOps
{
	// Non-zeros of consecutive rows of a matrix under construction
	struct RowChunk
	{
		vector<SparseMatrixIndex> RowEnds; // end of each row in Cols and Values
		vector<SparseMatrixIndex> Cols;
		vector<double> Values;

		void Add(SparseMatrixIndex col, double value)
		{
			Cols.push_back(col);
			Values.push_back(value);
		}

		void EndRow()
		{
			RowEnds.push_back(Cols.size());
		}
	};

	// Builds the nRows x nCols matrix whose rows are computed by computeRows(firstRow, lastRow, chunk), which must add to
	// the chunk the non-zeros of the rows firstRow, ..., lastRow - 1, in this order, sorted by column, and call
	// chunk.EndRow() after each row. The rows are split into chunks processed in parallel. The chunk boundaries are
	// multiples of rowBlockSize, so that computeRows can process whole blocks of rows.
	template <class ComputeRows>
	SparseMatrix BuildByRows(BigNumber nRows, BigNumber nCols, int rowBlockSize, ComputeRows computeRows)
	{
		assert(nRows % rowBlockSize == 0);
		BigNumber nBlocks = nRows / rowBlockSize;
		BigNumber nChunks = min(nBlocks, (BigNumber)(8 * Parallelism::NThreads())); // several chunks per thread, for load balancing
		vector<RowChunk> chunks(nChunks);
		auto chunkFirstRow = [&](BigNumber c) { return (nBlocks * c / nChunks) * rowBlockSize; };

		#pragma omp parallel for schedule(dynamic)
		for (BigNumber c = 0; c < nChunks; ++c)
		{
			chunks[c].RowEnds.reserve(chunkFirstRow(c + 1) - chunkFirstRow(c));
			computeRows(chunkFirstRow(c), chunkFirstRow(c + 1), chunks[c]);
			assert((BigNumber)chunks[c].RowEnds.size() == chunkFirstRow(c + 1) - chunkFirstRow(c));
		}

		// Concatenation of the chunks
		vector<SparseMatrixIndex> chunkOffsets(nChunks + 1, 0);
		for (BigNumber c = 0; c < nChunks; ++c)
			chunkOffsets[c + 1] = chunkOffsets[c] + chunks[c].Cols.size();

		SparseMatrix M(nRows, nCols);
		M.resizeNonZeros(chunkOffsets[nChunks]);
		#pragma omp parallel for schedule(dynamic)
		for (BigNumber c = 0; c < nChunks; ++c)
		{
			const RowChunk& chunk = chunks[c];
			BigNumber firstRow = chunkFirstRow(c);
			for (BigNumber r = 0; r < chunk.RowEnds.size(); ++r)
				M.outerIndexPtr()[firstRow + r + 1] = chunkOffsets[c] + chunk.RowEnds[r];
			copy(chunk.Cols.begin(), chunk.Cols.end(), M.innerIndexPtr() + chunkOffsets[c]);
			copy(chunk.Values.begin(), chunk.Values.end(), M.valuePtr() + chunkOffsets[c]);
		}
		M.outerIndexPtr()[0] = 0;
		return M;
	}

	// A * B, with the same result as Eigen's product: each coefficient is summed in the order of the non-zeros of the
	// row of A, the explicit zeros are kept, and the rows are sorted by column.
	inline SparseMatrix Multiply(const SparseMatrix& A, const SparseMatrix& B)
	{
		assert(A.cols() == B.rows());
		// For each thread: position of each column in the row being computed (-1 if absent), and that row
		ThreadLocal<vector<SparseMatrixIndex>> positions((size_t)B.cols(), (SparseMatrixIndex)-1);
		ThreadLocal<vector<pair<SparseMatrixIndex, double>>> rows;

		return BuildByRows(A.rows(), B.cols(), 1, [&](BigNumber firstRow, BigNumber lastRow, RowChunk& chunk)
		{
			vector<SparseMatrixIndex>& pos = positions.Local();
			vector<pair<SparseMatrixIndex, double>>& row = rows.Local();
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				// A single coefficient 1 in the row of A (e.g. an injection in a prolongation): the row of B, exactly
				SparseMatrix::InnerIterator a0(A, i);
				bool singleOne = false;
				if (a0 && a0.value() == 1.0)
				{
					SparseMatrix::InnerIterator a1 = a0;
					++a1;
					singleOne = !a1;
				}
				if (singleOne)
				{
					for (SparseMatrix::InnerIterator b(B, a0.col()); b; ++b)
						chunk.Add(b.col(), b.value() * 1.0);
					chunk.EndRow();
					continue;
				}
				row.clear();
				for (SparseMatrix::InnerIterator a(A, i); a; ++a)
				{
					for (SparseMatrix::InnerIterator b(B, a.col()); b; ++b)
					{
						SparseMatrixIndex& p = pos[b.col()];
						if (p < 0)
						{
							p = row.size();
							row.push_back({ b.col(), b.value() * a.value() });
						}
						else
							row[p].second += b.value() * a.value();
					}
				}
				sort(row.begin(), row.end(), [](const pair<SparseMatrixIndex, double>& x, const pair<SparseMatrixIndex, double>& y) { return x.first < y.first; });
				for (const pair<SparseMatrixIndex, double>& coeff : row)
				{
					chunk.Add(coeff.first, coeff.second);
					pos[coeff.first] = -1;
				}
				chunk.EndRow();
			}
		});
	}

	// A^T, with the same result as Eigen's: the coefficients of each row of A^T in the order of the rows of A, the explicit
	// zeros kept. In parallel: each thread counts, then copies, the coefficients of a contiguous range of rows of A.
	// A: a SparseMatrix (of any scalar type), or a map of row-major arrays (e.g. those of a column-major matrix, i.e. of
	// its transpose).
	template <class RowMajorSparseMatrix>
	SparseMatrixOf<typename RowMajorSparseMatrix::Scalar> Transpose(const RowMajorSparseMatrix& A)
	{
		static_assert(RowMajorSparseMatrix::IsRowMajor, "Transpose() iterates over the rows of A");
		using Result = SparseMatrixOf<typename RowMajorSparseMatrix::Scalar>;
		BigNumber nRows = A.rows(), nCols = A.cols();
		int nChunks = (int)min((BigNumber)Parallelism::NThreads(), max(nRows, (BigNumber)1));
		if (nChunks <= 1 || A.nonZeros() < 20000)
			return Result(A.transpose());
		// count[c][j]: number of coefficients of column j in the rows of chunk c, then the position of the first of them
		vector<vector<SparseMatrixIndex>> count(nChunks);
		#pragma omp parallel for schedule(static, 1)
		for (int c = 0; c < nChunks; ++c)
		{
			count[c].assign(nCols, 0);
			auto [begin, end] = Parallelism::Chunk(nRows, c, nChunks);
			for (BigNumber i = begin; i < end; ++i)
				for (typename RowMajorSparseMatrix::InnerIterator it(A, i); it; ++it)
					count[c][it.col()]++;
		}
		Result T(nCols, nRows);
		T.resizeNonZeros(A.nonZeros());
		SparseMatrixIndex* outer = T.outerIndexPtr();
		outer[0] = 0;
		for (BigNumber j = 0; j < nCols; ++j)
		{
			SparseMatrixIndex position = outer[j];
			for (int c = 0; c < nChunks; ++c)
			{
				SparseMatrixIndex n = count[c][j];
				count[c][j] = position;
				position += n;
			}
			outer[j + 1] = position;
		}
		#pragma omp parallel for schedule(static, 1)
		for (int c = 0; c < nChunks; ++c)
		{
			auto [begin, end] = Parallelism::Chunk(nRows, c, nChunks);
			for (BigNumber i = begin; i < end; ++i)
			{
				for (typename RowMajorSparseMatrix::InnerIterator it(A, i); it; ++it)
				{
					SparseMatrixIndex p = count[c][it.col()]++;
					T.innerIndexPtr()[p] = i;
					T.valuePtr()[p] = it.value();
				}
			}
		}
		return T;
	}

	// Symmetric matrix whose lower triangular part (diagonal included) is that of A: same matrix as
	// SparseMatrix(A.selfadjointView<Eigen::Lower>()), with the rows sorted by column.
	inline SparseMatrix FullFromLower(const SparseMatrix& A)
	{
		assert(A.rows() == A.cols());
		SparseMatrix At = Transpose(A); // row i of At: column i of A, whose part under the diagonal gives the part of row i above it
		return BuildByRows(A.rows(), A.cols(), 1, [&](BigNumber firstRow, BigNumber lastRow, RowChunk& chunk)
		{
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				for (SparseMatrix::InnerIterator it(A, i); it; ++it)
				{
					if (it.col() <= (SparseMatrixIndex)i)
						chunk.Add(it.col(), it.value());
				}
				for (SparseMatrix::InnerIterator it(At, i); it; ++it)
				{
					if (it.col() > (SparseMatrixIndex)i)
						chunk.Add(it.col(), it.value());
				}
				chunk.EndRow();
			}
		});
	}

	// Matrix whose i-th block of rows (blocks of blockSize rows) is that of A if fromA[i], that of B otherwise.
	// The coefficients of absolute value <= dropTolerance are dropped.
	inline SparseMatrix SelectRows(const vector<bool>& fromA, const SparseMatrix& A, const SparseMatrix& B, int blockSize, double dropTolerance)
	{
		assert(A.rows() == B.rows() && A.cols() == B.cols() && (BigNumber)fromA.size() * blockSize == (BigNumber)A.rows());
		return BuildByRows(A.rows(), A.cols(), blockSize, [&](BigNumber firstRow, BigNumber lastRow, RowChunk& chunk)
		{
			for (BigNumber i = firstRow; i < lastRow; ++i)
			{
				const SparseMatrix& M = fromA[i / blockSize] ? A : B;
				for (SparseMatrix::InnerIterator it(M, i); it; ++it)
				{
					if (abs(it.value()) > dropTolerance)
						chunk.Add(it.col(), it.value());
				}
				chunk.EndRow();
			}
		});
	}

	// A in precision U (e.g. the matrices of the levels of U-AMG in single precision): its structure copied, its values
	// rounded, in parallel (Eigen's cast() inserts the non-zeros one by one)
	template <typename U, typename T>
	SparseMatrixOf<U> Cast(const SparseMatrixOf<T>& A)
	{
		if (!A.isCompressed())
			return A.template cast<U>();
		SparseMatrixOf<U> M(A.rows(), A.cols());
		M.resizeNonZeros(A.nonZeros());
		copy(A.outerIndexPtr(), A.outerIndexPtr() + A.outerSize() + 1, M.outerIndexPtr());
		const SparseMatrixIndex* inner = A.innerIndexPtr();
		const T* values = A.valuePtr();
		SparseMatrixIndex* MInner = M.innerIndexPtr();
		U* MValues = M.valuePtr();
		BigNumber nnz = A.nonZeros();
		#pragma omp parallel for if (nnz > 20000)
		for (BigNumber p = 0; p < nnz; p++)
		{
			MInner[p] = inner[p];
			MValues[p] = (U)values[p];
		}
		return M;
	}
}
