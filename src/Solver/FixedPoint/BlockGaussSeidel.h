#pragma once
#include "BlockSOR.h"
using namespace std;

class BlockGaussSeidel : public IterativeSolver
{
protected:
	int _blockSize;
	Direction _direction;
	bool _hybrid = false;

	vector<Eigen::FullPivLU<DenseMatrix>> invD;
	RowMajorSparseMatrix _rowMajorA;
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
		if (!A.IsRowMajor)
			this->_rowMajorA = A;

		auto nb = A.rows() / _blockSize;
		this->invD = vector<Eigen::FullPivLU<DenseMatrix>>(nb);

		NumberParallelLoop<EmptyResultChunk> parallelLoop(nb);
		parallelLoop.Execute([this, &A](BigNumber i, ParallelChunk<EmptyResultChunk>* chunk)
			{
				DenseMatrix Di = A.block(i * _blockSize, i * _blockSize, _blockSize, _blockSize);
				this->invD[i].compute(Di);
			});

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
			{
				for (BigNumber i = 0; i < nb; ++i)
					ProcessBlockRow(i, b, x);
			}
			else if (_direction == Direction::Backward)
			{
				for (BigNumber i = 0; i < nb; ++i)
					ProcessBlockRow(nb - i - 1, b, x);
			}
			else if (_direction == Direction::Symmetric)
			{
				// Forward
				for (BigNumber i = 0; i < nb; ++i)
					ProcessBlockRow(i, b, x);
				// Backward
				for (BigNumber i = 0; i < nb; ++i)
					ProcessBlockRow(nb - i - 1, b, x);
			}
			else if (_direction == Direction::AlternatingForwardFirst || _direction == Direction::AlternatingBackwardFirst)
			{
				int modulo = _direction == Direction::AlternatingForwardFirst ? 0 : 1;
				if (this->IterationCount % 2 == modulo)
				{
					// Forward
					for (BigNumber i = 0; i < nb; ++i)
						ProcessBlockRow(i, b, x);
				}
				else
				{
					// Backward
					for (BigNumber i = 0; i < nb; ++i)
						ProcessBlockRow(nb - i - 1, b, x);
				}
			}
			else
				Utils::FatalError("direction not managed");
		}
		else // Hybrid
		{
			if (_direction == Direction::Forward)
			{
				NumberParallelLoop<EmptyResultChunk> parallelLoop(nb);
				parallelLoop.ExecuteChunk([this, &b, &x](ParallelChunk<EmptyResultChunk>* chunk)
					{
						for (BigNumber i = chunk->Start; i < chunk->End; i++)
							ProcessBlockRow(i, b, x);
					});
			}
			else if (_direction == Direction::Backward)
			{
				NumberParallelLoop<EmptyResultChunk> parallelLoop(nb);
				parallelLoop.ExecuteChunk([this, &b, &x, &nb](ParallelChunk<EmptyResultChunk>* chunk)
					{
						for (BigNumber i = chunk->Start; i < chunk->End; i++)
							ProcessBlockRow(nb - i - 1, b, x);
					});
			}
			else if (_direction == Direction::Symmetric)
			{
				// Forward
				NumberParallelLoop<EmptyResultChunk> parallelLoop(nb);
				parallelLoop.ExecuteChunk([this, &b, &x](ParallelChunk<EmptyResultChunk>* chunk)
					{
						for (BigNumber i = chunk->Start; i < chunk->End; i++)
							ProcessBlockRow(i, b, x);
					});
				// Backward
				parallelLoop.ExecuteChunk([this, &b, &x, &nb](ParallelChunk<EmptyResultChunk>* chunk)
					{
						for (BigNumber i = chunk->Start; i < chunk->End; i++)
							ProcessBlockRow(nb - i - 1, b, x);
					});
			}
			else if (_direction == Direction::AlternatingForwardFirst || _direction == Direction::AlternatingBackwardFirst)
			{
				NumberParallelLoop<EmptyResultChunk> parallelLoop(nb);
				int modulo = _direction == Direction::AlternatingForwardFirst ? 0 : 1;
				if (this->IterationCount % 2 == modulo)
				{
					// Forward
					parallelLoop.ExecuteChunk([this, &b, &x](ParallelChunk<EmptyResultChunk>* chunk)
						{
							for (BigNumber i = chunk->Start; i < chunk->End; i++)
								ProcessBlockRow(i, b, x);
						});
				}
				else
				{
					// Backward
					parallelLoop.ExecuteChunk([this, &b, &x, &nb](ParallelChunk<EmptyResultChunk>* chunk)
					{
						for (BigNumber i = chunk->Start; i < chunk->End; i++)
							ProcessBlockRow(nb - i - 1, b, x);
					});
				}
			}
			else
				Utils::FatalError("direction not managed");
		}

		result.SetX(x);
		double sweepWork = 2 * A.nonZeros() + nb * pow(_blockSize, 2);
		result.AddWorkInFlops(_direction == Direction::Symmetric ? 2*sweepWork : sweepWork);
		return result;
	}

	void ProcessBlockRow(BigNumber currentBlockRow, const Vector& b, Vector& x)
	{
		const SparseMatrix& A = *this->Matrix;

		// BlockRow i: [ --- Li --- | Di | --- Ui --- ]

		/*
		auto nb = A.rows() / _blockSize;
		// x_new                                               x_new                                  x_old                                               x_old
		x.segment(i * _blockSize, _blockSize) = -_omega * Li * x.head(i * _blockSize) - _omega * Ui * x.tail((nb - i - 1) * _blockSize) + (1 - _omega)*Di*x.segment(i * _blockSize, _blockSize) + _omega * bi;
		x.segment(i * _blockSize, _blockSize) = this->invD.block(i * _blockSize, 0, _blockSize, _blockSize) * x.segment(i * _blockSize, _blockSize);*/

		Vector tmp_x = b.segment(currentBlockRow * _blockSize, _blockSize);

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
				if (iBlock != jBlock) // Li and Ui
					tmp_x(k) -= a_ij * x(j);
			}
		}

		x.segment(currentBlockRow * _blockSize, _blockSize) = this->invD[currentBlockRow].solve(tmp_x);
	}
};