#pragma once
#include "IterativeSolver.h"
#include "../Utils/ParallelLoop.h"
using namespace std;

class BlockDampedJacobi : public IterativeSolver
{
protected:
	int _blockSize;
	double _omega;

	vector<Eigen::FullPivLU<Eigen::MatrixXd>> invD;
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
		this->invD = vector<Eigen::FullPivLU<Eigen::MatrixXd>>(nb);

		NumberParallelLoop<EmptyResultChunk> parallelLoop(nb);
		parallelLoop.Execute([this](BigNumber i, ParallelChunk<EmptyResultChunk>* chunk)
			{
				Eigen::MatrixXd Di = this->A.block(i * _blockSize, i * _blockSize, _blockSize, _blockSize);
				this->invD[i].compute(Di);
			});
	}

private:
	Eigen::VectorXd ExecuteOneIteration(const Eigen::VectorXd& b, Eigen::VectorXd& xOld) override
	{
		auto nb = A.rows() / _blockSize;

		Eigen::VectorXd xNew(xOld.rows());

		NumberParallelLoop<EmptyResultChunk> parallelLoop(nb);
		parallelLoop.Execute([this, b, xOld, &xNew](BigNumber i, ParallelChunk<EmptyResultChunk>* chunk)
			{
				ProcessBlockRow(i, b, xOld, xNew);
			});
			
		return xNew;
	}

	inline void ProcessBlockRow(BigNumber currentBlockRow, const Eigen::VectorXd& b, const Eigen::VectorXd& xOld, Eigen::VectorXd& xNew)
	{
		// BlockRow i: [ --- Li --- | Di | --- Ui --- ]

		Eigen::VectorXd tmp_x = _omega * b.segment(currentBlockRow * _blockSize, _blockSize);

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

	virtual void Serialize(ostream& os) const override
	{
		os << "Damped Jacobi (omega=" << _omega << ")";
	}
};

class BlockJacobi : public BlockDampedJacobi
{
public:
	BlockJacobi(int blockSize) : BlockDampedJacobi(blockSize, 1) {}

	virtual void Serialize(ostream& os) const override
	{
		os << "Block Jacobi (blockSize=" << _blockSize << ")";
	}
};

class Jacobi : public BlockJacobi
{
public: 
	Jacobi() : BlockJacobi(1) {}

	virtual void Serialize(ostream& os) const override
	{
		os << "Jacobi";
	}
};