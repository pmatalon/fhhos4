#pragma once
#include "BlockSOR.h"
#include "BlockDiagonalSolver.h"
#include "DiagonalBlockPositions.h"
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
public:
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

		auto nb = A.rows() / _blockSize;
		this->diagBlockSolver.Setup(A, _blockSize);
		_diagBlock.Setup(A, _blockSize);
		SetupSamePattern();

		this->SetupComputationalWork = nb * 2.0/3.0*pow(_blockSize, 3)*1e-6;
	}

private:
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
				HybridSweep(b, x, false);
			else if (_direction == Direction::Backward)
				HybridSweep(b, x, true);
			else if (_direction == Direction::Symmetric)
			{
				HybridSweep(b, x, false); // Forward
				HybridSweep(b, x, true);  // Backward
			}
			else if (_direction == Direction::AlternatingForwardFirst || _direction == Direction::AlternatingBackwardFirst)
			{
				int modulo = _direction == Direction::AlternatingForwardFirst ? 0 : 1;
				if (this->IterationCount % 2 == modulo)
					HybridSweep(b, x, false); // Forward
				else
					HybridSweep(b, x, true);  // Backward
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

	void Sweep(const Vector& b, Vector& x, bool backward, bool xEquals0)
	{
		BigNumber nb = this->Matrix->rows() / _blockSize;
		Vector tmp_x(_blockSize), tmp_xi(_blockSize); // work vectors for ProcessBlockRow()
		OffDiagonalBlocks blocks = !xEquals0 ? OffDiagonalBlocks::Both : (backward ? OffDiagonalBlocks::Upper : OffDiagonalBlocks::Lower);
		for (BigNumber i = 0; i < nb; ++i)
			ProcessBlockRow(backward ? nb - i - 1 : i, b, x, tmp_x, tmp_xi, blocks, backward);
	}

	// Each thread runs a Gauss-Seidel sweep on its own contiguous chunk of block rows.
	// The result depends on the chunks, i.e. on the number of threads (schedule(static)), and on the timing of the
	// threads: x is read by a thread while another one writes it (it varies from one run to the next).
	void HybridSweep(const Vector& b, Vector& x, bool backward)
	{
		BigNumber nb = this->Matrix->rows() / _blockSize;
		#pragma omp parallel
		{
			Vector tmp_x(_blockSize), tmp_xi(_blockSize); // work vectors, allocated once per thread
			#pragma omp for schedule(static)
			for (BigNumber i = 0; i < nb; i++)
				ProcessBlockRow(backward ? nb - i - 1 : i, b, x, tmp_x, tmp_xi, OffDiagonalBlocks::Both, backward);
		}
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

		// Solved in a work vector: solved directly in x, its intermediate values (e.g. Eigen zeroes the destination
		// of a product, then accumulates) would be read by the other threads in the hybrid version.
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