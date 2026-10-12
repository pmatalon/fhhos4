#pragma once
#include "../../Utils/Utils.h"
using namespace std;

// Chunks of the hybrid Gauss-Seidel sweeps (GaussSeidel and BlockGaussSeidel with hybrid = true). The sweep sequence of
// the block rows (their natural order, or SweepOrder) is split into contiguous chunks, one per thread. Each thread sweeps
// its chunk by Gauss-Seidel, and the chunks are coupled by Jacobi: a block row reads the current values of its own
// chunk, and those of the other chunks from the start of the sweep. The result only depends on the number of chunks,
// not on the timing of the threads; with one chunk, it is the sequential sweep.
// Only the boundary block rows, which have a non-zero in another chunk, read the other chunks. Before the sweep, the
// threads copy the values of x that they read (those of the exported block rows) into a vector of old values
// (SaveExported()); after a barrier, the sweep processes the boundary block rows with their external non-zeros times the
// old values and their internal ones times the current values (BoundaryRowRhs()). The other block rows, the interior
// ones, keep the kernels of the sequential sweeps.
// If the chunks are poorly separated (a numbering without locality), most rows are boundary rows, and the Jacobi
// between the chunks makes the sweep a poor smoother, even a divergent one: when the external non-zeros are more than
// MaxExternalFraction of the non-zeros, the sweeps are the sequential ones (Parallel() false).
// The non-zeros of each row of a boundary block row are described by ranges of positions in the row: the external
// ones, then the internal ones before the diagonal block, and after it. In the natural order, a chunk is a range of
// columns (the rows are sorted by column): at most two external ranges and one internal range on each side.
// Set up once per matrix, number of chunks and order: with an order, one pass over the matrix; in the natural order, the
// first and last non-zeros of the rows, and the boundary rows.
class HybridChunks
{
public:
	// Which internal products a boundary row subtracts: from x = 0, those by the values not computed yet are skipped
	// (they are 0), as in the sequential sweeps: the upper ones in a forward sweep, the lower ones in a backward sweep
	enum class InternalProducts { Both, Lower, Upper };

private:
	using Range = pair<SparseMatrixIndex, SparseMatrixIndex>; // positions [first, second), from the start of the row

	struct Chunk
	{
		BigNumber Begin = 0, End = 0; // positions [Begin, End) in the sweep sequence
		vector<BigNumber> Boundary;   // positions of the boundary block rows, increasing
		// Row k of the boundary block row q (r = q * blockSize + k), in Ranges: external ranges [RangeBegin[r],
		// LowerBegin[r]), internal ones before the diagonal block [LowerBegin[r], UpperBegin[r]), after it
		// [UpperBegin[r], RangeBegin[r + 1])
		vector<BigNumber> RangeBegin;
		vector<BigNumber> LowerBegin;
		vector<BigNumber> UpperBegin;
		vector<Range> Ranges;
		BigNumber NExternal = 0;      // number of external non-zeros
	};

	const SparseMatrix* _A = nullptr;
	int _blockSize = 1;
	const vector<BigNumber>* _order = nullptr;
	int _nChunks = 0; // 0: not set up
	bool _parallel = false;
	vector<Chunk> _chunks;
	vector<BigNumber> _exported; // block rows read by another chunk, increasing
	Vector _old;                 // the values of x at the start of the sweep, on the exported block rows

public:
	// Above this fraction of external non-zeros, the sweeps are sequential. Measured (U-AMG, FCG or GCR, 4 threads): with
	// 32% on the fine level of fhhos4's Cube-tet (GMSH numbering), the hybrid sweeps diverge; with 15% on its level 2,
	// and 3-7% on HArDCore3D's cubes (4-8 threads), they converge as well as the sequential ones (PERFORMANCE.md).
	static constexpr double MaxExternalFraction = 0.2;

	// Whether the sweeps run in parallel: on several chunks, well separated
	bool Parallel() const
	{
		return _parallel;
	}

	int Count() const
	{
		return (int)_chunks.size();
	}

	bool IsSetUp(int nChunks, const vector<BigNumber>* order) const
	{
		return _nChunks == nChunks && _order == order;
	}

	// To call when the matrix changes
	void Clear()
	{
		_nChunks = 0;
		_parallel = false;
		_chunks.clear();
		_exported.clear();
		_old = Vector();
	}

	// Block row at a position of the sweep sequence
	BigNumber BlockRow(BigNumber position) const
	{
		return _order ? (*_order)[position] : position;
	}

	void Setup(const SparseMatrix& A, int blockSize, const vector<BigNumber>* order, int nChunks)
	{
		Clear();
		_A = &A;
		_blockSize = blockSize;
		_order = order;
		_nChunks = nChunks;
		if (nChunks == 1)
			return;
		_chunks.assign(nChunks, Chunk());
		BigNumber nb = A.rows() / blockSize;
		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();

		for (int c = 0; c < nChunks; c++)
			tie(_chunks[c].Begin, _chunks[c].End) = Parallelism::Chunk(nb, c, nChunks);
		// The chunk of each block row, with an order. In the natural order, a chunk is a range of block columns: a block
		// row is interior if the first and last columns of each of its rows are in it.
		vector<int> chunkOf(order ? nb : 0);
		if (order)
		{
			#pragma omp parallel for
			for (int c = 0; c < nChunks; c++)
			{
				for (BigNumber k = _chunks[c].Begin; k < _chunks[c].End; k++)
					chunkOf[BlockRow(k)] = c;
			}
		}

		#pragma omp parallel for schedule(dynamic)
		for (int c = 0; c < nChunks; c++)
		{
			Chunk& chunk = _chunks[c];
			auto inChunk = [&](BigNumber J) { return order ? chunkOf[J] == c : J >= chunk.Begin && J < chunk.End; };
			chunk.RangeBegin.push_back(0);
			for (BigNumber k = chunk.Begin; k < chunk.End; k++)
			{
				BigNumber I = BlockRow(k);
				bool boundary = false;
				if (order)
				{
					for (SparseMatrixIndex p = outer[I * blockSize]; p < outer[(I + 1) * blockSize] && !boundary; p++)
						boundary = !inChunk(col[p] / blockSize);
				}
				else
				{
					for (BigNumber i = I * blockSize; i < (I + 1) * blockSize && !boundary; i++)
						boundary = outer[i] < outer[i + 1] && (!inChunk(col[outer[i]] / blockSize) || !inChunk(col[outer[i + 1] - 1] / blockSize));
				}
				if (!boundary)
					continue;
				chunk.Boundary.push_back(k);
				for (BigNumber i = I * blockSize; i < (I + 1) * blockSize; i++)
				{
					// Ranges of the non-zeros of row i selected by isSelected(J), J their block column
					auto addRanges = [&](auto isSelected)
					{
						bool open = false; // whether the last range can be extended
						for (SparseMatrixIndex p = outer[i]; p < outer[i + 1]; p++)
						{
							bool selected = isSelected((BigNumber)col[p] / blockSize);
							if (selected && open)
								chunk.Ranges.back().second++;
							else if (selected)
								chunk.Ranges.push_back({ p - outer[i], p - outer[i] + 1 });
							open = selected;
						}
					};
					addRanges([&](BigNumber J) { return !inChunk(J); });
					chunk.LowerBegin.push_back(chunk.Ranges.size());
					addRanges([&](BigNumber J) { return J < I && inChunk(J); });
					chunk.UpperBegin.push_back(chunk.Ranges.size());
					addRanges([&](BigNumber J) { return J > I && inChunk(J); });
					chunk.RangeBegin.push_back(chunk.Ranges.size());
					for (BigNumber j = chunk.RangeBegin[chunk.RangeBegin.size() - 2]; j < chunk.LowerBegin.back(); j++)
						chunk.NExternal += chunk.Ranges[j].second - chunk.Ranges[j].first;
				}
			}
		}

		BigNumber nExternal = 0;
		for (const Chunk& chunk : _chunks)
			nExternal += chunk.NExternal;
		if (nExternal > MaxExternalFraction * A.nonZeros())
		{
			_chunks.clear();
			return;
		}
		_parallel = true;

		// The exported block rows: the block columns of the external non-zeros
		vector<char> exported(nb, false);
		for (const Chunk& chunk : _chunks)
		{
			BigNumber nBoundary = chunk.Boundary.size();
			for (BigNumber q = 0; q < nBoundary; q++)
			{
				for (BigNumber i = BlockRow(chunk.Boundary[q]) * blockSize, r = q * blockSize; r < (q + 1) * blockSize; i++, r++)
				{
					for (BigNumber j = chunk.RangeBegin[r]; j < chunk.LowerBegin[r]; j++)
					{
						for (SparseMatrixIndex p = outer[i] + chunk.Ranges[j].first; p < outer[i] + chunk.Ranges[j].second; p++)
							exported[col[p] / blockSize] = true;
					}
				}
			}
		}
		_exported.clear();
		for (BigNumber J = 0; J < nb; J++)
		{
			if (exported[J])
				_exported.push_back(J);
		}
		_old.resize(A.rows());
	}

	// In a parallel region: f(c) for the chunks of the calling thread (its own one if the region has one thread per chunk)
	template <class F>
	void ForEachChunk(F f) const
	{
		for (int c = Parallelism::ThreadNumber(); c < Count(); c += Parallelism::TeamSize())
			f(c);
	}

	// Before the sweep, by every thread of the parallel region (then a barrier before any thread writes x): copies the
	// values of x of the exported block rows, split evenly among the threads, into the old values
	void SaveExported(const Vector& x)
	{
		auto [begin, end] = Parallelism::ThreadChunk(_exported.size());
		for (BigNumber e = begin; e < end; e++)
		{
			BigNumber i0 = _exported[e] * _blockSize;
			for (BigNumber i = i0; i < i0 + _blockSize; i++)
				_old[i] = x[i];
		}
	}

	// Sweeps the chunk c: interior(begin, end) for the runs of interior block rows (positions [begin, end) of the sweep
	// sequence, processed backward if backward), boundary(q, blockRow) for its boundary block rows, in the sweep order
	template <class Interior, class Boundary>
	void SweepChunk(int c, bool backward, Interior interior, Boundary boundary) const
	{
		const Chunk& chunk = _chunks[c];
		BigNumber nBoundary = chunk.Boundary.size();
		if (!backward)
		{
			BigNumber position = chunk.Begin;
			for (BigNumber q = 0; q < nBoundary; q++)
			{
				interior(position, chunk.Boundary[q]);
				boundary(q, BlockRow(chunk.Boundary[q]));
				position = chunk.Boundary[q] + 1;
			}
			interior(position, chunk.End);
		}
		else
		{
			BigNumber position = chunk.End;
			for (BigNumber q = nBoundary; q-- > 0; )
			{
				interior(chunk.Boundary[q] + 1, position);
				boundary(q, BlockRow(chunk.Boundary[q]));
				position = chunk.Boundary[q];
			}
			interior(chunk.Begin, position);
		}
	}

	// b_i - sum_j a_ij x_j over the j outside the diagonal block, for the row i = blockRow * blockSize + k of the boundary
	// block row q of the chunk c: the external products with the old values (SaveExported(); none if xEquals0: x = 0),
	// the internal ones selected by 'products' with the current values of x
	double BoundaryRowRhs(int c, BigNumber q, BigNumber blockRow, int k, const Vector& b, const Vector& x, bool xEquals0, InternalProducts products) const
	{
		const Chunk& chunk = _chunks[c];
		const SparseMatrixIndex* col = _A->innerIndexPtr();
		const double* val = _A->valuePtr();
		const double* xp = x.data();
		const double* old = _old.data();
		BigNumber i = blockRow * _blockSize + k;
		BigNumber r = q * _blockSize + k;
		SparseMatrixIndex start = _A->outerIndexPtr()[i];
		double s = b[i];
		if (!xEquals0)
		{
			for (BigNumber j = chunk.RangeBegin[r]; j < chunk.LowerBegin[r]; j++)
			{
				for (SparseMatrixIndex p = start + chunk.Ranges[j].first; p < start + chunk.Ranges[j].second; p++)
					s -= val[p] * old[col[p]];
			}
		}
		BigNumber first = products == InternalProducts::Upper ? chunk.UpperBegin[r] : chunk.LowerBegin[r];
		BigNumber last = products == InternalProducts::Lower ? chunk.UpperBegin[r] : chunk.RangeBegin[r + 1];
		for (BigNumber j = first; j < last; j++)
		{
			for (SparseMatrixIndex p = start + chunk.Ranges[j].first; p < start + chunk.Ranges[j].second; p++)
				s -= val[p] * xp[col[p]];
		}
		return s;
	}

	// The same for all the rows of the boundary block row: rhs[k] for its row k. samePattern: its rows have the same
	// non-zero pattern: processed together, as in BlockGaussSeidel::BlockRowRhsSamePattern().
	void BoundaryRowRhs(int c, BigNumber q, BigNumber blockRow, const Vector& b, const Vector& x, bool xEquals0, InternalProducts products, double* rhs, bool samePattern) const
	{
		if (samePattern)
		{
			switch (_blockSize)
			{
				case 2: BoundaryRowRhsSamePattern<2>(c, q, blockRow, b, x, xEquals0, products, rhs); return;
				case 3: BoundaryRowRhsSamePattern<3>(c, q, blockRow, b, x, xEquals0, products, rhs); return;
				case 4: BoundaryRowRhsSamePattern<4>(c, q, blockRow, b, x, xEquals0, products, rhs); return;
				case 5: BoundaryRowRhsSamePattern<5>(c, q, blockRow, b, x, xEquals0, products, rhs); return;
				case 6: BoundaryRowRhsSamePattern<6>(c, q, blockRow, b, x, xEquals0, products, rhs); return;
			}
		}
		for (int k = 0; k < _blockSize; k++)
			rhs[k] = BoundaryRowRhs(c, q, blockRow, k, b, x, xEquals0, products);
	}

private:
	// The rows of a block row with the same pattern are stored one after the other with the same length, and have the
	// same ranges (those of its first row): processed together, column by column
	template <int BlockSize>
	void BoundaryRowRhsSamePattern(int c, BigNumber q, BigNumber blockRow, const Vector& b, const Vector& x, bool xEquals0, InternalProducts products, double* rhs) const
	{
		const Chunk& chunk = _chunks[c];
		BigNumber i0 = blockRow * BlockSize;
		BigNumber r0 = q * BlockSize;
		SparseMatrixIndex start = _A->outerIndexPtr()[i0];
		SparseMatrixIndex length = _A->outerIndexPtr()[i0 + 1] - start;
		const SparseMatrixIndex* col = _A->innerIndexPtr() + start;
		const double* val[BlockSize]; // coefficients of each row
		double s[BlockSize];
		for (int k = 0; k < BlockSize; k++)
		{
			val[k] = _A->valuePtr() + start + k * length;
			s[k] = b[i0 + k];
		}
		if (!xEquals0)
		{
			for (BigNumber j = chunk.RangeBegin[r0]; j < chunk.LowerBegin[r0]; j++)
			{
				for (SparseMatrixIndex p = chunk.Ranges[j].first; p < chunk.Ranges[j].second; p++)
				{
					double old_j = _old[col[p]];
					for (int k = 0; k < BlockSize; k++)
						s[k] -= val[k][p] * old_j;
				}
			}
		}
		const double* xp = x.data();
		BigNumber first = products == InternalProducts::Upper ? chunk.UpperBegin[r0] : chunk.LowerBegin[r0];
		BigNumber last = products == InternalProducts::Lower ? chunk.UpperBegin[r0] : chunk.RangeBegin[r0 + 1];
		for (BigNumber j = first; j < last; j++)
		{
			for (SparseMatrixIndex p = chunk.Ranges[j].first; p < chunk.Ranges[j].second; p++)
			{
				double x_j = xp[col[p]];
				for (int k = 0; k < BlockSize; k++)
					s[k] -= val[k][p] * x_j;
			}
		}
		for (int k = 0; k < BlockSize; k++)
			rhs[k] = s[k];
	}
};
