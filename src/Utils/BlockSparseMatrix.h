#pragma once
#include <algorithm>
#include <type_traits>
#include <vector>
#include "Types.h"
#include "Parallelism.h"
#include "Utils.h"
using namespace std;

// Calls f(std::integral_constant<int, B>()) with B = blockSize if the kernels are compiled for this size (the faces of
// degree 1 to 4: 2, 3, 4, 5 in 2D, 3, 6, 10, 15 in 3D, and of degree 1 to 3 in Q^k in 3D: 4, 9, 16), with B = 0
// otherwise (the sizes read at run time, any of them, also different from a block to another). With the size known at
// compile time, the products of the small blocks are unrolled and vectorized: the kernels of block CSR are 1.5-2.5x
// faster than in CSR, without it no faster (PERFORMANCE.md).
template <class F>
void DispatchBlockSize(int blockSize, F f)
{
	switch (blockSize)
	{
		case 2:  f(integral_constant<int, 2>());  return;
		case 3:  f(integral_constant<int, 3>());  return;
		case 4:  f(integral_constant<int, 4>());  return;
		case 5:  f(integral_constant<int, 5>());  return;
		case 6:  f(integral_constant<int, 6>());  return;
		case 9:  f(integral_constant<int, 9>());  return;
		case 10: f(integral_constant<int, 10>()); return;
		case 15: f(integral_constant<int, 15>()); return;
		case 16: f(integral_constant<int, 16>()); return;
		default: f(integral_constant<int, 0>());  return;
	}
}

// A square sparse matrix stored by dense blocks (block CSR): its rows and columns are grouped in blocks of consecutive
// DoFs (the DoFs of the faces of a hybrid discretization), each pair of coupled blocks stored as a dense block,
// column-major, with one column index per block instead of one per coefficient: per coefficient, 4 + 4/b^2 bytes in
// float and 8 + 4/b^2 in double for blocks of size b, instead of 8 and 12 in CSR, which the products, bound by the
// memory bandwidth, read faster. The blocks have the sizes given by BlockStart (variable, e.g. variable degrees), all
// BlockSize when they are the same. The coefficients of a block missing in the CSR matrix are stored as zeros. The
// products sum each row in the order of its columns, with one accumulator, as the CSR kernels of the smoothers: the same
// results, bit for bit (Eigen 5's sparse product sums with two accumulators: the same to rounding).
// In precision T: double, or float (the levels of U-AMG in mixed precision), the sums in double.
template <typename T>
class BlockSparseMatrix
{
public:
	vector<SparseMatrixIndex> BlockStart; // the first row (and column) of each block, then the number of rows
	int BlockSize = 0;                    // the size of all the blocks, 0 if they have different sizes
	vector<SparseMatrixIndex> RowStart;   // the blocks of the block row I: RowStart[I], ..., RowStart[I + 1] - 1
	vector<SparseMatrixIndex> Cols;       // the block column of each block, increasing in each block row
	vector<size_t> ValueStart;            // the values of the block row I, from ValueStart[I]: its blocks one after the other
	vector<T> Values;

	// The largest block (the kernels sum the rows of a block row in an array of this size on the stack)
	static constexpr int MaxBlockSize = 256;

	BigNumber NBlocks() const { return BlockStart.empty() ? 0 : BlockStart.size() - 1; }
	BigNumber Rows() const { return BlockStart.empty() ? 0 : BlockStart.back(); }
	size_t NonZeros() const { return Values.size(); } // the coefficients stored
	int Size(BigNumber I) const { return BlockStart[I + 1] - BlockStart[I]; }
	int MaxSize() const
	{
		int max = BlockSize;
		for (BigNumber I = 0; I < NBlocks() && BlockSize == 0; I++)
			max = std::max(max, Size(I));
		return max;
	}

	// From the CSR matrix A (of any precision) whose rows and columns are grouped in blocks of blockSize DoFs
	template <typename S>
	static BlockSparseMatrix FromCSR(const SparseMatrixOf<S>& A, int blockSize)
	{
		vector<SparseMatrixIndex> blockStart(A.rows() / blockSize + 1);
		for (size_t I = 0; I < blockStart.size(); I++)
			blockStart[I] = I * blockSize;
		return FromCSR(A, std::move(blockStart));
	}

	// From the CSR matrix A (of any precision) whose rows and columns are grouped in the blocks starting at blockStart.
	// In parallel, by block rows.
	template <typename S>
	static BlockSparseMatrix FromCSR(const SparseMatrixOf<S>& A, vector<SparseMatrixIndex> blockStart)
	{
		assert(A.rows() == A.cols() && blockStart.back() == A.rows());
		BlockSparseMatrix M;
		M.BlockStart = std::move(blockStart);
		BigNumber nb = M.NBlocks();
		M.BlockSize = nb > 0 ? M.Size(0) : 0;
		for (BigNumber I = 1; I < nb && M.BlockSize > 0; I++)
			if (M.Size(I) != M.BlockSize)
				M.BlockSize = 0;
		if (M.MaxSize() > MaxBlockSize)
			Utils::FatalError("The block CSR storage supports blocks of at most " + to_string(MaxBlockSize) + " DoFs.");
		// The block of each column (variable sizes)
		vector<SparseMatrixIndex> blockOf;
		if (M.BlockSize == 0)
		{
			blockOf.resize(A.cols());
			#pragma omp parallel for
			for (BigNumber J = 0; J < nb; J++)
				for (SparseMatrixIndex j = M.BlockStart[J]; j < M.BlockStart[J + 1]; j++)
					blockOf[j] = J;
		}
		auto block = [&](SparseMatrixIndex col) { return M.BlockSize > 0 ? col / M.BlockSize : blockOf[col]; };

		// The block columns of the block row I: the union of those of its rows, increasing (usually the same in all of
		// them: the first one's, checked)
		auto blockCols = [&](BigNumber I, vector<SparseMatrixIndex>& cols, vector<SparseMatrixIndex>& row, vector<SparseMatrixIndex>& merged)
		{
			cols.clear();
			for (SparseMatrixIndex i = M.BlockStart[I]; i < M.BlockStart[I + 1]; i++)
			{
				row.clear();
				for (typename SparseMatrixOf<S>::InnerIterator it(A, i); it; ++it)
				{
					SparseMatrixIndex J = block(it.col());
					if (row.empty() || row.back() != J)
						row.push_back(J);
				}
				if (i == M.BlockStart[I])
					cols.swap(row);
				else if (row != cols)
				{
					merged.clear();
					set_union(cols.begin(), cols.end(), row.begin(), row.end(), back_inserter(merged));
					cols.swap(merged);
				}
			}
		};

		vector<SparseMatrixIndex> nBlocks(nb);
		vector<size_t> nValues(nb);
		#pragma omp parallel
		{
			vector<SparseMatrixIndex> cols, row, merged;
			#pragma omp for schedule(dynamic, 1024)
			for (BigNumber I = 0; I < nb; I++)
			{
				blockCols(I, cols, row, merged);
				nBlocks[I] = cols.size();
				size_t width = 0;
				for (SparseMatrixIndex J : cols)
					width += M.Size(J);
				nValues[I] = (size_t)M.Size(I) * width;
			}
		}
		M.RowStart.assign(nb + 1, 0);
		M.ValueStart.assign(nb + 1, 0);
		for (BigNumber I = 0; I < nb; I++)
		{
			M.RowStart[I + 1] = M.RowStart[I] + nBlocks[I];
			M.ValueStart[I + 1] = M.ValueStart[I] + nValues[I];
		}
		M.Cols.resize(M.RowStart[nb]);
		M.Values.resize(M.ValueStart[nb]);
		#pragma omp parallel
		{
			vector<SparseMatrixIndex> cols, row, merged;
			vector<size_t> offsets; // of the blocks of the block row in its values
			#pragma omp for schedule(dynamic, 1024)
			for (BigNumber I = 0; I < nb; I++)
			{
				blockCols(I, cols, row, merged);
				copy(cols.begin(), cols.end(), M.Cols.begin() + M.RowStart[I]);
				int n = M.Size(I);
				T* values = M.Values.data() + M.ValueStart[I];
				fill(values, values + nValues[I], T(0));
				offsets.resize(cols.size());
				for (size_t p = 0, offset = 0; p < cols.size(); offset += (size_t)n * M.Size(cols[p]), p++)
					offsets[p] = offset;
				for (int r = 0; r < n; r++)
				{
					size_t p = 0; // the block of the current coefficient: both increasing
					for (typename SparseMatrixOf<S>::InnerIterator it(A, M.BlockStart[I] + r); it; ++it)
					{
						SparseMatrixIndex J = block(it.col());
						while (cols[p] != J)
							p++;
						values[offsets[p] + (size_t)(it.col() - M.BlockStart[J]) * n + r] = (T)it.value();
					}
				}
			}
		}
		return M;
	}

	// acc = rIn_I - (A v)_I on the block row I (rIn null: - (A v)_I), summed in double in the order of the columns. B: the
	// size of the blocks if it is known at compile time (all the blocks), 0 otherwise (read at run time). acc: Size(I)
	// values.
	template <int B>
	inline void RowResidual(BigNumber I, const T* v, const T* rIn, double* acc) const
	{
		const int n = B > 0 ? B : Size(I);
		SparseMatrixIndex first = B > 0 ? I * B : BlockStart[I];
		for (int r = 0; r < n; r++)
			acc[r] = rIn ? (double)rIn[first + r] : 0.0;
		const T* block = Values.data() + ValueStart[I];
		for (SparseMatrixIndex p = RowStart[I]; p < RowStart[I + 1]; p++)
		{
			SparseMatrixIndex J = Cols[p];
			const int m = B > 0 ? B : Size(J);
			const T* vJ = v + (B > 0 ? J * B : BlockStart[J]);
			for (int c = 0; c < m; c++)
			{
				double vc = vJ[c];
				const T* column = block + c * n;
				for (int r = 0; r < n; r++)
					acc[r] -= (double)column[r] * vc;
			}
			block += n * m;
		}
	}

	// y = A x
	void Multiply(const T* x, T* y) const
	{
		DispatchBlockSize(BlockSize, [&](auto size)
		{
			constexpr int B = decltype(size)::value;
			#pragma omp parallel if (NonZeros() > 20000)
			{
				double acc[B > 0 ? B : MaxBlockSize];
				#pragma omp for schedule(static)
				for (BigNumber I = 0; I < NBlocks(); I++)
				{
					RowResidual<B>(I, x, nullptr, acc);
					SparseMatrixIndex first = B > 0 ? I * B : BlockStart[I];
					for (int r = 0; r < (B > 0 ? B : Size(I)); r++)
						y[first + r] = (T)-acc[r];
				}
			}
		});
	}

	VectorOf<T> operator*(const VectorOf<T>& x) const
	{
		VectorOf<T> y(Rows());
		Multiply(x.data(), y.data());
		return y;
	}

	// The diagonal block of the block row I (dense, in double), or zero if it is not stored
	DenseMatrix DiagonalBlock(BigNumber I) const
	{
		int n = Size(I);
		DenseMatrix D = DenseMatrix::Zero(n, n);
		const T* block = Values.data() + ValueStart[I];
		for (SparseMatrixIndex p = RowStart[I]; p < RowStart[I + 1]; p++)
		{
			if (Cols[p] == I)
			{
				for (int c = 0; c < n; c++)
					for (int r = 0; r < n; r++)
						D(r, c) = block[c * n + r];
				break;
			}
			block += (size_t)n * Size(Cols[p]);
		}
		return D;
	}
};
