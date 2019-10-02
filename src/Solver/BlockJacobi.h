#pragma once
#include "IterativeSolver.h"
#include "../Utils/ParallelLoop.h"
using namespace std;

class BlockDampedJacobi : public IterativeSolver
{
protected:
	int _blockSize;
	double _omega;

	vector<Eigen::FullPivLU<DenseMatrix>> invD;
	Eigen::SparseMatrix<double, Eigen::RowMajor> _rowMajorA;
public:

	BlockDampedJacobi(int blockSize, double omega)
	{
		this->_blockSize = blockSize;
		this->_omega = omega;
	}

	virtual void Serialize(ostream& os) const override
	{
		os << "Block-Damped Jacobi (blockSize=" << _blockSize << ", omega=" << _omega << ")";
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
		parallelLoop.Execute([this](BigNumber i, ParallelChunk<EmptyResultChunk>* chunk)
			{
				DenseMatrix Di = this->A.block(i * _blockSize, i * _blockSize, _blockSize, _blockSize);
				this->invD[i].compute(Di);
			});

		this->SetupComputationalWork = 2.0 / 3.0*pow(_blockSize, 3);
	}

private:
	IterationResult ExecuteOneIteration(const Vector& b, Vector& xOld, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);

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

	inline void ProcessBlockRow(BigNumber currentBlockRow, const Vector& b, const Vector& xOld, Vector& xNew)
	{
		// BlockRow i: [ --- Li --- | Di | --- Ui --- ]

		Vector tmp_x = _omega * b.segment(currentBlockRow * _blockSize, _blockSize);

		for (int k = 0; k < _blockSize; k++)
		{
			BigNumber iBlock = currentBlockRow;
			BigNumber i = iBlock * _blockSize + k;
			// RowMajor --> the following line iterates over the non-zeros of the i-th row.
			for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A.IsRowMajor ? A : _rowMajorA, i); it; ++it)
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

class DampedJacobi : public BlockDampedJacobi
{
public: 
	DampedJacobi(double omega) : BlockDampedJacobi(1, omega) {}

	static string Code() { return "dj"; };

	virtual void Serialize(ostream& os) const override
	{
		os << "Damped Jacobi (omega=" << _omega << ")";
	}
};

class BlockJacobi : public BlockDampedJacobi
{
public:
	BlockJacobi(int blockSize) : BlockDampedJacobi(blockSize, 1) {}

	static string Code() { return "bj"; };

	virtual void Serialize(ostream& os) const override
	{
		os << "Block Jacobi (blockSize=" << _blockSize << ")";
	}
};

class Jacobi : public BlockJacobi
{
public: 
	Jacobi() : BlockJacobi(1) {}

	static string Code() { return "j"; };

	virtual void Serialize(ostream& os) const override
	{
		os << "Jacobi";
	}
};