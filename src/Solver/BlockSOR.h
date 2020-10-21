#pragma once
#include "IterativeSolver.h"
#include "EigenSparseLU.h"
using namespace std;

enum class Direction : unsigned
{
	Forward = 0,
	Backward = 1
};

class BlockSOR : public IterativeSolver
{
protected:
	int _blockSize;
	double _omega;
	Direction _direction;

	vector<Eigen::FullPivLU<DenseMatrix>> invD;
	RowMajorSparseMatrix _rowMajorA;
public:
	BlockSOR(int blockSize, double omega) : BlockSOR(blockSize, omega, Direction::Forward) {}

	BlockSOR(int blockSize, double omega, Direction direction)
	{
		this->_blockSize = blockSize;
		this->_omega = omega;
		this->_direction = direction;
	}

	virtual void Serialize(ostream& os) const override
	{
		if (_blockSize != 1)
			os << "Block ";
		if (_omega == 1)
			os << "Gauss-Seidel ";
		else
			os << "SOR ";

		os << " (";
		if (_blockSize != 1)
			os << "blockSize=" << _blockSize << ", ";
		if (_omega != 1)
			os << "omega=" << _omega << ", ";
		os << "direction=" << (_direction == Direction::Forward ? "forward" : "backward") << ")";
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

		this->SetupComputationalWork = 2.0 / 3.0*pow(_blockSize, 3);
	}

private:
	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);

		const SparseMatrix& A = *this->Matrix;

		auto nb = A.rows() / _blockSize;
		if (_direction == Direction::Forward)
		{
			for (BigNumber i = 0; i < nb; ++i)
				ProcessBlockRow(i, b, x);
		}
		else
		{
			for (BigNumber i = 0; i < nb; ++i)
				ProcessBlockRow(nb - i - 1, b, x);
		}

		result.SetX(x);
		result.AddCost(2 * A.nonZeros() + nb * pow(_blockSize, 2));
		return result;
	}

	inline void ProcessBlockRow(BigNumber currentBlockRow, const Vector& b, Vector& x)
	{
		const SparseMatrix& A = *this->Matrix;

		// BlockRow i: [ --- Li --- | Di | --- Ui --- ]

		/*
		auto nb = A.rows() / _blockSize;
		// x_new                                               x_new                                  x_old                                               x_old
		x.segment(i * _blockSize, _blockSize) = -_omega * Li * x.head(i * _blockSize) - _omega * Ui * x.tail((nb - i - 1) * _blockSize) + (1 - _omega)*Di*x.segment(i * _blockSize, _blockSize) + _omega * bi;
		x.segment(i * _blockSize, _blockSize) = this->invD.block(i * _blockSize, 0, _blockSize, _blockSize) * x.segment(i * _blockSize, _blockSize);*/

		Vector tmp_x = _omega * b.segment(currentBlockRow * _blockSize, _blockSize);

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
					tmp_x(k) += (1 - _omega) * a_ij * x(j);
				else // Li and Ui
					tmp_x(k) += -_omega * a_ij * x(j);
			}
		}

		x.segment(currentBlockRow * _blockSize, _blockSize) = this->invD[currentBlockRow].solve(tmp_x);
	}
};