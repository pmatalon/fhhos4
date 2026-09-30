#pragma once
#include "../../Utils/Types.h"
#include "../../Utils/ParallelLoop.h"
using namespace std;

enum class DiagBlockSolveMethod : unsigned
{
	// The LU factorizations (with full pivoting) of the blocks are computed in the setup: each solve is a real (triangular) solve.
	LU,
	// The inverses of the blocks are computed in the setup: each solve is then a small matrix-vector product.
	// Slightly faster (up to ~20% on the whole smoothing), but less accurate if the blocks are ill-conditioned.
	Inverse
};

// Solves the linear systems D_i x_i = b_i, where the D_i are the diagonal blocks of a matrix.
// Used in the block smoothers, where it is done for each block row at each iteration: Solve() does not allocate memory.
class BlockDiagonalSolver
{
public:
	inline static DiagBlockSolveMethod DefaultMethod = DiagBlockSolveMethod::LU;

private:
	int _blockSize = 0;
	DiagBlockSolveMethod _method = DiagBlockSolveMethod::LU;
	// blockSize x (nBlocks * blockSize): the i-th block is the inverse of D_i (Inverse),
	// or its LU factors, as stored by Eigen::FullPivLU (LU)
	DenseMatrix _blocks;
	// LU: row and column permutations (P and Q such that D_i = P^-1 L U Q^-1), and rank of each block
	vector<int> _permP;
	vector<int> _permQ;
	vector<int> _rank;
public:
	void Setup(const SparseMatrix& A, int blockSize, DiagBlockSolveMethod method = DefaultMethod)
	{
		_blockSize = blockSize;
		_method = method;
		BigNumber nb = A.rows() / blockSize;
		_blocks = DenseMatrix(blockSize, nb * blockSize);
		if (_method == DiagBlockSolveMethod::LU)
		{
			_permP = vector<int>(nb * blockSize);
			_permQ = vector<int>(nb * blockSize);
			_rank = vector<int>(nb);
		}

		NumberParallelLoop<EmptyResultChunk> parallelLoop(nb);
		parallelLoop.Execute([this, &A](BigNumber i)
			{
				DenseMatrix Di = A.block(i * _blockSize, i * _blockSize, _blockSize, _blockSize);
				Eigen::FullPivLU<DenseMatrix> lu(Di);
				auto block = _blocks.middleCols(i * _blockSize, _blockSize);
				if (_method == DiagBlockSolveMethod::Inverse)
					block = lu.inverse();
				else
				{
					block = lu.matrixLU();
					for (int k = 0; k < _blockSize; k++)
					{
						_permP[i * _blockSize + k] = lu.permutationP().indices()(k);
						_permQ[i * _blockSize + k] = lu.permutationQ().indices()(k);
					}
					_rank[i] = lu.nonzeroPivots();
				}
			});
	}

	// x = D_i^{-1} b (b can have several columns).
	// b is used as work space: its content is destroyed. x must not alias b.
	template <class RhsT, class DestT>
	void Solve(BigNumber i, Eigen::MatrixBase<RhsT>& b, Eigen::MatrixBase<DestT>& x) const
	{
		auto Di = _blocks.middleCols(i * _blockSize, _blockSize);
		if (_method == DiagBlockSolveMethod::Inverse)
			x.derived().noalias() = Di * b;
		else
		{
			// Same steps as Eigen::FullPivLU::solve(), without its temporaries
			const int* P = &_permP[i * _blockSize];
			const int* Q = &_permQ[i * _blockSize];
			int rank = _rank[i];
			// x = P b
			for (int k = 0; k < _blockSize; k++)
				x.row(P[k]) = b.row(k);
			// x = L^-1 x
			Di.template triangularView<Eigen::UnitLower>().solveInPlace(x);
			// x = U^-1 x (non-singular part)
			Di.topLeftCorner(rank, rank).template triangularView<Eigen::Upper>().solveInPlace(x.topRows(rank));
			// x = Q x (via b)
			for (int k = 0; k < rank; k++)
				b.row(Q[k]) = x.row(k);
			for (int k = rank; k < _blockSize; k++)
				b.row(Q[k]).setZero();
			x = b;
		}
	}
};
