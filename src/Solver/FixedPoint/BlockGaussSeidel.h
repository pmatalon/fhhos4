#pragma once
#include "BlockSOR.h"
#include "BlockDiagonalSolver.h"
#include "DiagonalBlockPositions.h"
#include "HybridChunks.h"
using namespace std;

class BlockGaussSeidel : public IterativeSolver
{
protected:
	int _blockSize;
	Direction _direction;
	bool _hybrid = false;

	BlockDiagonalSolver diagBlockSolver; // solves the systems with the diagonal blocks
	DiagonalBlockPositions _diagBlock;   // positions of the diagonal block in each row
	// For each block row, whether its rows have the same non-zero pattern (always the case for the HHO matrices)
	vector<char> _samePattern;
	HybridChunks _chunks; // of the hybrid sweeps (see Chunks())
public:
	// If set, every sweep (forward or backward, pre- or post-smoothing) processes the block rows in this order: the
	// downwind order of a convection-dominated matrix (DownwindOrder()), which all the sweeps then follow
	const vector<BigNumber>* SweepOrder = nullptr;

	BlockGaussSeidel(int blockSize, Direction direction, bool hybrid)
	{
		this->_blockSize = blockSize;
		this->_direction = direction;
		this->_hybrid = hybrid;
	}

	virtual void Serialize(ostream& os) const override
	{
		if (_hybrid)
			os << "Hybrid ";
		if (_direction == Direction::Symmetric)
			os << "Symmetric ";
		os << "Block Gauss-Seidel";

		if (_direction != Direction::Symmetric || _blockSize != 1)
			os << " (";
		if (_blockSize != 1)
		{
			os << "blockSize=" << _blockSize;
			if (_direction != Direction::Symmetric)
				os << ", ";
		}
		if (_direction == Direction::Forward)
			os << "direction=forward";
		else if (_direction == Direction::Backward)
			os << "direction=backward";
		else if (_direction == Direction::AlternatingForwardFirst)
			os << "alternating directions, forward first";
		else if (_direction == Direction::AlternatingBackwardFirst)
			os << "alternating directions, backward first";

		if (_direction != Direction::Symmetric || _blockSize != 1)
			os << ")";
	}

	//------------------------------------------------//
	// Big assumption: the matrix must be symmetric!  //
	//------------------------------------------------//

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		SetupLocal(A);
	}

	// On several processes: Gauss-Seidel in the rows of the process, Jacobi between the processes (see SweepRhs())
	void Setup(const DistributedMatrix& A) override
	{
		SetDistributedMatrix(A);
		SetupLocal(A.Diag);
	}

private:
	void SetupLocal(const SparseMatrix& A)
	{
		auto nb = A.rows() / _blockSize;
		this->diagBlockSolver.Setup(A, _blockSize);
		_diagBlock.Setup(A, _blockSize);
		SetupSamePattern();
		_chunks.Clear();
		if (_hybrid)
			Chunks();

		this->SetupComputationalWork = nb * 2.0/3.0*pow(_blockSize, 3)*1e-6;
	}

public:
	// The factorized diagonal blocks, the positions of the blocks and, with the same order, the hybrid chunks of
	// another block Gauss-Seidel set up for A with the same block size
	bool SetupFrom(const SparseMatrix& A, const IterativeSolver& other) override
	{
		const BlockGaussSeidel* o = dynamic_cast<const BlockGaussSeidel*>(&other);
		if (!o || o->Matrix != &A || o->Distributed || o->_blockSize != _blockSize || this->ComputeExactSolution)
			return false;
		IterativeSolver::Setup(A);
		CopySetup(*o);
		return true;
	}

	bool SetupFrom(const DistributedMatrix& A, const IterativeSolver& other) override
	{
		const BlockGaussSeidel* o = dynamic_cast<const BlockGaussSeidel*>(&other);
		if (!o || o->Distributed != &A || o->_blockSize != _blockSize || this->ComputeExactSolution)
			return false;
		SetDistributedMatrix(A);
		CopySetup(*o);
		return true;
	}

private:
	void CopySetup(const BlockGaussSeidel* o)
	{
		diagBlockSolver = o->diagBlockSolver;
		_diagBlock = o->_diagBlock;
		_samePattern = o->_samePattern;
		_chunks.Clear();
		if (_hybrid)
		{
			if (o->_hybrid && o->SweepOrder == SweepOrder)
				_chunks = o->_chunks;
			else
				Chunks();
		}
		this->SetupComputationalWork = o->SetupComputationalWork;
	}

	void CopySetup(const BlockGaussSeidel& o)
	{
		CopySetup(&o);
	}

	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);
		assert(!computeResidual && !computeAx);

		const SparseMatrix& A = *this->Matrix;

		auto nb = A.rows() / _blockSize;

		if (!_hybrid)
		{
			if (_direction == Direction::Forward)
				Sweep(b, x, false, xEquals0);
			else if (_direction == Direction::Backward)
				Sweep(b, x, true, xEquals0);
			else if (_direction == Direction::Symmetric)
			{
				Sweep(b, x, false, xEquals0); // Forward
				Sweep(b, x, true, false);     // Backward
			}
			else if (_direction == Direction::AlternatingForwardFirst || _direction == Direction::AlternatingBackwardFirst)
			{
				int modulo = _direction == Direction::AlternatingForwardFirst ? 0 : 1;
				if (this->IterationCount % 2 == modulo)
					Sweep(b, x, false, xEquals0); // Forward
				else
					Sweep(b, x, true, xEquals0);  // Backward
			}
			else
				Utils::FatalError("direction not managed");
		}
		else // Hybrid
		{
			if (_direction == Direction::Forward)
				HybridSweep(b, x, false, xEquals0);
			else if (_direction == Direction::Backward)
				HybridSweep(b, x, true, xEquals0);
			else if (_direction == Direction::Symmetric)
			{
				HybridSweep(b, x, false, xEquals0); // Forward
				HybridSweep(b, x, true, false);     // Backward
			}
			else if (_direction == Direction::AlternatingForwardFirst || _direction == Direction::AlternatingBackwardFirst)
			{
				int modulo = _direction == Direction::AlternatingForwardFirst ? 0 : 1;
				if (this->IterationCount % 2 == modulo)
					HybridSweep(b, x, false, xEquals0); // Forward
				else
					HybridSweep(b, x, true, xEquals0);  // Backward
			}
			else
				Utils::FatalError("direction not managed");
		}

		xEquals0 = false;

		result.SetX(x);
		double sweepWork = 2 * A.nonZeros() + nb * pow(_blockSize, 2);
		result.AddWorkInFlops(_direction == Direction::Symmetric ? 2*sweepWork : sweepWork);
		return result;
	}

	// Which off-diagonal blocks of a block row are multiplied by x: in the first sweep from x = 0, those of the block
	// rows not processed yet multiply zeros (the upper ones in a forward sweep, the lower ones in a backward sweep)
	enum class OffDiagonalBlocks { Both, Lower, Upper };

	void Sweep(const Vector& givenB, Vector& x, bool backward, bool xEquals0)
	{
		Vector distributedRhs;
		const Vector& b = SweepRhs(givenB, x, xEquals0, distributedRhs);
		BigNumber nb = this->Matrix->rows() / _blockSize;
		Vector tmp_x(_blockSize), tmp_xi(_blockSize); // work vectors for ProcessBlockRow()
		if (SweepOrder && xEquals0)
			x = Vector::Zero(this->Matrix->rows());
		SweepRange(b, x, backward, xEquals0, 0, nb, tmp_x, tmp_xi);
	}

	// The block rows at the positions [begin, end) of the sweep sequence: the natural order (processed backward if
	// backward), or SweepOrder (all the sweeps follow it)
	void SweepRange(const Vector& b, Vector& x, bool backward, bool xEquals0, BigNumber begin, BigNumber end, Vector& tmp_x, Vector& tmp_xi)
	{
		if (SweepOrder)
		{
			// The block rows not processed yet are not the upper ones: from x = 0, all the off-diagonal blocks are
			// multiplied (by zeros for those)
			for (BigNumber k = begin; k < end; ++k)
				ProcessBlockRow((*SweepOrder)[k], b, x, tmp_x, tmp_xi, OffDiagonalBlocks::Both, true);
			return;
		}
		OffDiagonalBlocks blocks = !xEquals0 ? OffDiagonalBlocks::Both : (backward ? OffDiagonalBlocks::Upper : OffDiagonalBlocks::Lower);
		if (backward)
		{
			for (BigNumber i = end; i-- > begin; )
				ProcessBlockRow(i, b, x, tmp_x, tmp_xi, blocks, true);
		}
		else
		{
			for (BigNumber i = begin; i < end; ++i)
				ProcessBlockRow(i, b, x, tmp_x, tmp_xi, blocks, false);
		}
	}

	// Hybrid sweep (see HybridChunks): each thread sweeps its own contiguous chunk of the block rows (or of SweepOrder) by
	// Gauss-Seidel, and reads the values of the other chunks from the start of the sweep (Jacobi between the chunks). The
	// interior block rows are processed as in Sweep(). The result depends on the number of threads, not on their
	// timing; on one thread, or if the chunks are poorly separated (HybridChunks::Parallel()), it is Sweep(). From x = 0,
	// the products by the values of the other chunks and by those not computed yet are skipped, as in Sweep() (except in
	// SweepOrder).
	void HybridSweep(const Vector& givenB, Vector& x, bool backward, bool xEquals0)
	{
		if (!Chunks().Parallel())
		{
			Sweep(givenB, x, backward, xEquals0);
			return;
		}
		Vector distributedRhs;
		const Vector& b = SweepRhs(givenB, x, xEquals0, distributedRhs);
		BigNumber n = this->Matrix->rows();
		if (xEquals0 && (SweepOrder || x.size() != n))
			x = Vector::Zero(n); // as Sweep()
		bool interleaveRows = backward || SweepOrder; // as in Sweep()
		bool backwardChunks = backward && !SweepOrder;
		auto products = !xEquals0 || SweepOrder ? HybridChunks::InternalProducts::Both : (backward ? HybridChunks::InternalProducts::Upper : HybridChunks::InternalProducts::Lower);
		#pragma omp parallel num_threads(_chunks.Count())
		{
			Vector tmp_x(_blockSize), tmp_xi(_blockSize); // work vectors, allocated once per thread
			if (!xEquals0)
			{
				_chunks.SaveExported(x);
				#pragma omp barrier
			}
			_chunks.ForEachChunk([&](int c)
			{
				_chunks.SweepChunk(c, backwardChunks,
					[&](BigNumber begin, BigNumber end) { SweepRange(b, x, backward, xEquals0, begin, end, tmp_x, tmp_xi); },
					[&](BigNumber q, BigNumber blockRow)
					{
						_chunks.BoundaryRowRhs(c, q, blockRow, b, x, xEquals0, products, tmp_x.data(), interleaveRows && _samePattern[blockRow]);
						this->diagBlockSolver.Solve(blockRow, tmp_x, tmp_xi);
						x.segment(blockRow * _blockSize, _blockSize) = tmp_xi;
					});
			});
		}
	}

	// The chunks of the hybrid sweeps: one per thread. Set up with the matrix, again if the number of threads or
	// SweepOrder changes.
	const HybridChunks& Chunks()
	{
		int nChunks = Parallelism::NThreads();
		if (!_chunks.IsSetUp(nChunks, SweepOrder))
			_chunks.Setup(*this->Matrix, _blockSize, SweepOrder, nChunks);
		return _chunks;
	}

	void SetupSamePattern()
	{
		const SparseMatrix& A = *this->Matrix;
		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();
		BigNumber nb = A.rows() / _blockSize;
		_samePattern.resize(nb);
		#pragma omp parallel for
		for (BigNumber iBlock = 0; iBlock < nb; iBlock++)
		{
			BigNumber i0 = iBlock * _blockSize;
			SparseMatrixIndex length = outer[i0 + 1] - outer[i0];
			bool same = true;
			for (int k = 1; k < _blockSize && same; k++)
			{
				BigNumber i = i0 + k;
				same = outer[i + 1] - outer[i] == length && equal(col + outer[i], col + outer[i + 1], col + outer[i0]);
			}
			_samePattern[iBlock] = same;
		}
	}

	// interleaveRows: whether to use BlockRowRhsSamePattern() when possible. Measured on the HHO matrices (U-AMG, k=2):
	// in a forward sweep, BlockRowRhs() is faster (it reads the matrix as one sequential stream, which the hardware
	// prefetcher follows); in a backward sweep, BlockRowRhsSamePattern() is (BlockRowRhs() jumps back at each block row,
	// and the prefetcher doesn't follow: reading the rows together starts their memory accesses together).
	void ProcessBlockRow(BigNumber currentBlockRow, const Vector& b, Vector& x, Vector& tmp_x, Vector& tmp_xi, OffDiagonalBlocks blocks, bool interleaveRows)
	{
		// BlockRow i: [ --- Li --- | Di | --- Ui --- ]

		/*
		auto nb = A.rows() / _blockSize;
		// x_new                                               x_new                                  x_old                                               x_old
		x.segment(i * _blockSize, _blockSize) = -_omega * Li * x.head(i * _blockSize) - _omega * Ui * x.tail((nb - i - 1) * _blockSize) + (1 - _omega)*Di*x.segment(i * _blockSize, _blockSize) + _omega * bi;
		x.segment(i * _blockSize, _blockSize) = this->invD.block(i * _blockSize, 0, _blockSize, _blockSize) * x.segment(i * _blockSize, _blockSize);*/

		// tmp_x = b_i - (Li | Ui) * x
		if (interleaveRows && _samePattern[currentBlockRow])
		{
			switch (_blockSize)
			{
				case 2: BlockRowRhsSamePattern<2>(currentBlockRow, b, x, tmp_x, blocks); break;
				case 3: BlockRowRhsSamePattern<3>(currentBlockRow, b, x, tmp_x, blocks); break;
				case 4: BlockRowRhsSamePattern<4>(currentBlockRow, b, x, tmp_x, blocks); break;
				case 5: BlockRowRhsSamePattern<5>(currentBlockRow, b, x, tmp_x, blocks); break;
				case 6: BlockRowRhsSamePattern<6>(currentBlockRow, b, x, tmp_x, blocks); break;
				default: BlockRowRhs(currentBlockRow, b, x, tmp_x, blocks);
			}
		}
		else
			BlockRowRhs(currentBlockRow, b, x, tmp_x, blocks);

		this->diagBlockSolver.Solve(currentBlockRow, tmp_x, tmp_xi);
		x.segment(currentBlockRow * _blockSize, _blockSize) = tmp_xi;
	}

	// rhs = b_i - (Li | Ui) * x for the block row i: for each row, b_k - a_kj*x_j - ..., the products subtracted one by one
	// in the order of the row. The sum of a row is a chain of dependent operations: it is kept in a register.
	// The blocks excluded by 'blocks' multiply zeros: skipped (s - 0 = s, exactly, up to the sign of a zero s).
	void BlockRowRhs(BigNumber iBlock, const Vector& b, const Vector& x, Vector& rhs, OffDiagonalBlocks blocks) const
	{
		const SparseMatrix& A = *this->Matrix;
		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();
		const double* val = A.valuePtr();
		const double* xp = x.data();
		for (int k = 0; k < _blockSize; k++)
		{
			BigNumber i = iBlock * _blockSize + k;
			double s = b[i];
			if (blocks != OffDiagonalBlocks::Upper) // Li
			{
				for (SparseMatrixIndex p = outer[i]; p < _diagBlock.Begin[i]; p++)
					s -= val[p] * xp[col[p]];
			}
			if (blocks != OffDiagonalBlocks::Lower) // Ui
			{
				for (SparseMatrixIndex p = _diagBlock.End[i]; p < outer[i + 1]; p++)
					s -= val[p] * xp[col[p]];
			}
			rhs[k] = s;
		}
	}

	// Same as BlockRowRhs (same operations in the same order), when the rows of the block row have the same pattern:
	// they are then stored one after the other with the same length, and processed together, column by column.
	// Their sums are independent chains of operations, which the CPU interleaves.
	template <int BlockSize>
	void BlockRowRhsSamePattern(BigNumber iBlock, const Vector& b, const Vector& x, Vector& rhs, OffDiagonalBlocks blocks) const
	{
		const SparseMatrix& A = *this->Matrix;
		BigNumber i0 = iBlock * BlockSize;
		SparseMatrixIndex start = A.outerIndexPtr()[i0];
		SparseMatrixIndex length = A.outerIndexPtr()[i0 + 1] - start;
		SparseMatrixIndex diagBegin = _diagBlock.Begin[i0] - start;
		SparseMatrixIndex diagEnd = _diagBlock.End[i0] - start;
		const SparseMatrixIndex* col = A.innerIndexPtr() + start;
		const double* val[BlockSize]; // coefficients of each row
		double s[BlockSize];
		for (int k = 0; k < BlockSize; k++)
		{
			val[k] = A.valuePtr() + start + k * length;
			s[k] = b[i0 + k];
		}
		const double* xp = x.data();

		if (blocks != OffDiagonalBlocks::Upper) // Li
		{
			for (SparseMatrixIndex p = 0; p < diagBegin; p++)
			{
				double x_j = xp[col[p]];
				for (int k = 0; k < BlockSize; k++)
					s[k] -= val[k][p] * x_j;
			}
		}
		if (blocks != OffDiagonalBlocks::Lower) // Ui
		{
			for (SparseMatrixIndex p = diagEnd; p < length; p++)
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