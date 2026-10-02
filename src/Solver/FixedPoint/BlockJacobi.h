#pragma once
#include "../IterativeSolver.h"
#include "../../Utils/SparseMatrixOps.h"
#include "BlockDiagonalSolver.h"
using namespace std;

class BlockJacobi : public IterativeSolver
{
protected:
	int _blockSize;
	double _omega;

	BlockDiagonalSolver diagBlockSolver; // solves the systems with the diagonal blocks
	RowMajorSparseMatrix _rowMajorA;
public:

	BlockJacobi(int blockSize, double omega = 1)
	{
		assert(omega > 0 && omega < 2);
		this->_blockSize = blockSize;
		this->_omega = omega;
	}

	virtual void Serialize(ostream& os) const override
	{
		if (_blockSize == 1)
		{
			os << "Jacobi";
			if (_omega != 1)
			{
				os << " (omega=";
				if (_omega == 2.0/3.0)
					os << "2/3";
				else
					os << _omega;
				os << ")";
			}
		}
		else
		{
			os << "Block Jacobi (blockSize=" << _blockSize;
			if (_omega != 1)
			{
				os << ", omega=";
				if (_omega == 2.0 / 3.0)
					os << "2/3";
				else
				os <<  _omega;
			}
			os << ")";
		}
	}

	//------------------------------------------------//
	// Big assumption: the matrix must be symmetric!  //
	//------------------------------------------------//

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		if (!A.IsRowMajor)
			this->_rowMajorA = A;

		auto nb = A.rows() / _blockSize;
		this->diagBlockSolver.Setup(A, _blockSize);

		this->SetupComputationalWork = nb * 2.0/3.0*pow(_blockSize, 3)*1e-6;
	}

private:
	IterationResult ExecuteOneIteration(const Vector& b, Vector& xOld, bool& xEquals0, bool computeResidual, bool computeAx, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);
		assert(!computeResidual && !computeAx);

		const SparseMatrix& A = *this->Matrix;

		auto nb = A.rows() / _blockSize;

		Vector xNew(xOld.rows());

		#pragma omp parallel
		{
			Vector tmp_x(_blockSize); // work vector, allocated once per thread
			#pragma omp for
			for (BigNumber i = 0; i < nb; i++)
				ProcessBlockRow(i, b, xOld, xNew, tmp_x);
		}

		xOld.swap(xNew);
		result.SetX(xOld);
		result.AddWorkInFlops(2 * A.nonZeros() + nb * pow(_blockSize, 2));
		return result;
	}

protected:
	inline void ProcessBlockRow(BigNumber currentBlockRow, const Vector& b, const Vector& xOld, Vector& xNew, Vector& tmp_x)
	{
		const SparseMatrix& A = *this->Matrix;

		// BlockRow i: [ --- Li --- | Di | --- Ui --- ]

		tmp_x = _omega * b.segment(currentBlockRow * _blockSize, _blockSize);

		for (int k = 0; k < _blockSize; k++)
		{
			BigNumber iBlock = currentBlockRow;
			BigNumber i = iBlock * _blockSize + k;
			// RowMajor --> the following line iterates over the non-zeros of the i-th row.
			for (RowMajorSparseMatrix::InnerIterator it(A.IsRowMajor ? A : _rowMajorA, i); it; ++it)
			{
				auto j = it.col();
				auto jBlock = j / this->_blockSize;
				auto a_ij = it.value();
				if (iBlock == jBlock) // Di
					tmp_x(k) += (1 - _omega) * a_ij * xOld(j);
				else // Li and Ui
					tmp_x(k) += -_omega * a_ij * xOld(j);
			}
		}

		auto xNew_i = xNew.segment(currentBlockRow * _blockSize, _blockSize);
		this->diagBlockSolver.Solve(currentBlockRow, tmp_x, xNew_i);
	}

public:
	// J = I - omega * D^-1 * A
	SparseMatrix IterationMatrix()
	{
		return IterationMatrix(vector<bool>(this->Matrix->rows() / _blockSize, true));
	}

	// Rows of J of the block rows i such that blockRows[i] (the other rows are empty).
	// The coefficients of absolute value <= NonZeroCoefficients::ZeroThreshold are dropped.
	SparseMatrix IterationMatrix(const vector<bool>& blockRows)
	{
		const SparseMatrix& A = *this->Matrix;
		int bs = _blockSize;
		assert((BigNumber)blockRows.size() * bs == (BigNumber)A.rows());

		DenseMatrix oneMinusOmegaIdentity = (1 - _omega)*DenseMatrix::Identity(bs, bs);

		return SparseMatrixOps::BuildByRows(A.rows(), A.cols(), bs, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
		{
			vector<BigNumber> jBlocks; // block columns of the block row
			DenseMatrix blockRow;      // A_i, then J_i, stored densely by blocks
			DenseMatrix A_ij(bs, bs);
			DenseMatrix jacobiBlock(bs, bs);

			for (BigNumber iBlock = firstRow / bs; iBlock < lastRow / bs; ++iBlock)
			{
				if (!blockRows[iBlock])
				{
					for (int k = 0; k < bs; k++)
						chunk.EndRow();
					continue;
				}

				// Block row A_i: [ --- Li --- | Di | --- Ui --- ]
				jBlocks.clear();
				for (int k = 0; k < bs; k++)
				{
					for (RowMajorSparseMatrix::InnerIterator it(A, iBlock * bs + k); it; ++it)
						jBlocks.push_back(it.col() / bs);
				}
				sort(jBlocks.begin(), jBlocks.end());
				jBlocks.erase(unique(jBlocks.begin(), jBlocks.end()), jBlocks.end());

				blockRow.setZero(bs, jBlocks.size() * bs);
				for (int k = 0; k < bs; k++)
				{
					for (RowMajorSparseMatrix::InnerIterator it(A, iBlock * bs + k); it; ++it)
					{
						BigNumber b = lower_bound(jBlocks.begin(), jBlocks.end(), (BigNumber)it.col() / bs) - jBlocks.begin();
						blockRow(k, b * bs + it.col() % bs) = it.value();
					}
				}

				// J_i = [ --- -omega*Di^-1*Li --- | (1-omega)*I | --- -omega*Di^-1*Ui --- ]
				for (BigNumber b = 0; b < jBlocks.size(); b++)
				{
					if (jBlocks[b] == iBlock) // Di
						blockRow.middleCols(b * bs, bs) = oneMinusOmegaIdentity;
					else
					{
						A_ij = blockRow.middleCols(b * bs, bs);
						this->diagBlockSolver.Solve(iBlock, A_ij, jacobiBlock);
						jacobiBlock *= -_omega;
						blockRow.middleCols(b * bs, bs) = jacobiBlock;
					}
				}

				for (int k = 0; k < bs; k++)
				{
					for (BigNumber c = 0; c < blockRow.cols(); c++)
					{
						if (abs(blockRow(k, c)) > NonZeroCoefficients::ZeroThreshold)
							chunk.Add(jBlocks[c / bs] * bs + c % bs, blockRow(k, c));
					}
					chunk.EndRow();
				}
			}
		});
	}
};