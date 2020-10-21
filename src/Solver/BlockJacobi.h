#pragma once
#include "IterativeSolver.h"
#include "../Utils/ParallelLoop.h"
using namespace std;

class BlockJacobi : public IterativeSolver
{
protected:
	int _blockSize;
	double _omega;

	vector<Eigen::FullPivLU<DenseMatrix>> invD;
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
	IterationResult ExecuteOneIteration(const Vector& b, Vector& xOld, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);

		const SparseMatrix& A = *this->Matrix;

		auto nb = A.rows() / _blockSize;

		Vector xNew(xOld.rows());

		NumberParallelLoop<EmptyResultChunk> parallelLoop(nb);
		parallelLoop.Execute([this, b, xOld, &xNew](BigNumber i, ParallelChunk<EmptyResultChunk>* chunk)
			{
				ProcessBlockRow(i, b, xOld, xNew);
			});

		result.SetX(xNew);
		result.AddCost(2 * A.nonZeros() + nb * pow(_blockSize, 2));
		return result;
	}

protected:
	inline void ProcessBlockRow(BigNumber currentBlockRow, const Vector& b, const Vector& xOld, Vector& xNew)
	{
		const SparseMatrix& A = *this->Matrix;

		// BlockRow i: [ --- Li --- | Di | --- Ui --- ]

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
					tmp_x(k) += (1 - _omega) * a_ij * xOld(j);
				else // Li and Ui
					tmp_x(k) += -_omega * a_ij * xOld(j);
			}
		}

		xNew.segment(currentBlockRow * _blockSize, _blockSize) = this->invD[currentBlockRow].solve(tmp_x);
	}
};