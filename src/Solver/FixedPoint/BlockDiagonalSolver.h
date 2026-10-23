#pragma once
#include "../../Utils/Types.h"
#include "../../Utils/Parallelism.h"
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
	// or its LU factors, as stored by Eigen::FullPivLU (LU). In double, or in single precision (a matrix in single
	// precision: the solves of float vectors).
	DenseMatrix _blocks;
	Eigen::MatrixXf _singleBlocks;
	// LU: row and column permutations (P and Q such that D_i = P^-1 L U Q^-1), and rank of each block
	vector<int> _permP;
	vector<int> _permQ;
	vector<int> _rank;
public:
	// With 'blocks', only the diagonal blocks i such that (*blocks)[i] are factorized (the others can't be solved).
	// A in precision T: the factorizations in double, stored in T (float: for the solves of float vectors).
	template <typename T>
	void Setup(const SparseMatrixOf<T>& A, int blockSize, DiagBlockSolveMethod method = DefaultMethod, const vector<bool>* blocks = nullptr)
	{
		_blockSize = blockSize;
		_method = method;
		BigNumber nb = A.rows() / blockSize;
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& blocksT = BlocksIn<T>();
		blocksT.resize(blockSize, nb * blockSize);
		if (_method == DiagBlockSolveMethod::LU)
		{
			_permP = vector<int>(nb * blockSize);
			_permQ = vector<int>(nb * blockSize);
			_rank = vector<int>(nb);
		}

		#pragma omp parallel for
		for (BigNumber i = 0; i < nb; ++i)
		{
			if (blocks && !(*blocks)[i])
				continue;
			Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> DiT = A.block(i * _blockSize, i * _blockSize, _blockSize, _blockSize);
			DenseMatrix Di = DiT.template cast<double>();
			Eigen::FullPivLU<DenseMatrix> lu(Di);
			auto block = blocksT.middleCols(i * _blockSize, _blockSize);
			if (_method == DiagBlockSolveMethod::Inverse)
				block = DenseMatrix(lu.inverse()).template cast<T>();
			else
			{
				block = lu.matrixLU().template cast<T>();
				for (int k = 0; k < _blockSize; k++)
				{
					_permP[i * _blockSize + k] = lu.permutationP().indices()(k);
					_permQ[i * _blockSize + k] = lu.permutationQ().indices()(k);
				}
				_rank[i] = lu.nonzeroPivots();
			}
		}
	}

	// x = D_i^{-1} b (b can have several columns), in the precision of b and x (that of the matrix of the setup).
	// b is used as work space: its content is destroyed. x must not alias b.
	template <class RhsT, class DestT>
	void Solve(BigNumber i, Eigen::MatrixBase<RhsT>& b, Eigen::MatrixBase<DestT>& x) const
	{
		if constexpr (std::is_same_v<typename RhsT::Scalar, float>)
			Solve(_singleBlocks, i, b, x);
		else
			Solve(_blocks, i, b, x);
	}

private:
	template <typename T>
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& BlocksIn()
	{
		if constexpr (std::is_same_v<T, float>)
			return _singleBlocks;
		else
			return _blocks;
	}

	template <class Blocks, class RhsT, class DestT>
	void Solve(const Blocks& blocks, BigNumber i, Eigen::MatrixBase<RhsT>& b, Eigen::MatrixBase<DestT>& x) const
	{
		auto Di = blocks.middleCols(i * _blockSize, _blockSize);
		if (_method == DiagBlockSolveMethod::Inverse)
			x.derived().noalias() = Di * b;
		else
			SolveLU(Di, &_permP[i * _blockSize], &_permQ[i * _blockSize], _rank[i], b, x);
	}

public:

	// x = D^{-1} b, from the LU factors of D stored as by Eigen::FullPivLU (matrixLU(), the indices of permutationP()
	// and permutationQ(), nonzeroPivots()). b is used as work space: its content is destroyed. x must not alias b.
	template <class LU, class RhsT, class DestT>
	static void SolveLU(const LU& Di, const int* P, const int* Q, int rank, Eigen::MatrixBase<RhsT>& b, Eigen::MatrixBase<DestT>& x)
	{
		int blockSize = Di.rows();
		// Same steps as Eigen::FullPivLU::solve(), without its temporaries
		// x = P b
		for (int k = 0; k < blockSize; k++)
			x.row(P[k]) = b.row(k);
		// x = L^-1 x
		Di.template triangularView<Eigen::UnitLower>().solveInPlace(x);
		// x = U^-1 x (non-singular part)
		Di.topLeftCorner(rank, rank).template triangularView<Eigen::Upper>().solveInPlace(x.topRows(rank));
		// x = Q x (via b)
		for (int k = 0; k < rank; k++)
			b.row(Q[k]) = x.row(k);
		for (int k = rank; k < blockSize; k++)
			b.row(Q[k]).setZero();
		x = b;
	}
};
