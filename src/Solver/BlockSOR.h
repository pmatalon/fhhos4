#pragma once
#include <Eigen/Sparse>
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
	Eigen::SparseMatrix<double, Eigen::RowMajor> _rowMajorA;
public:

	BlockSOR(int blockSize, double omega, Direction direction)
	{
		this->_blockSize = blockSize;
		this->_omega = omega;
		this->_direction = direction;
	}

	virtual void Serialize(ostream& os) const override
	{
		os << "Block-SOR (blockSize=" << _blockSize << ", omega=" << _omega << ", direction=" << (_direction == Direction::Forward ? "forward" : "backward") << ")";
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
	IterationResult ExecuteOneIteration(const Vector& b, Vector& x, const IterationResult& oldResult) override
	{
		IterationResult result(oldResult);

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
			for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A.IsRowMajor ? A : _rowMajorA, i); it; ++it)
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

class SOR : public BlockSOR
{
public: 
	SOR(double omega) : BlockSOR(1, omega, Direction::Forward) {}

	virtual void Serialize(ostream& os) const override
	{
		os << "SOR (omega=" << _omega << ")";
	}
};

class ReverseSOR : public BlockSOR
{
public: 
	ReverseSOR(double omega) : BlockSOR(1, omega, Direction::Backward) {}

	virtual void Serialize(ostream& os) const override
	{
		os << "Reverse SOR (omega=" << _omega << ")";
	}
};

class GaussSeidel : public SOR
{
public: 
	GaussSeidel() : SOR(1) {}

	static string Code() { return "gs"; };

	virtual void Serialize(ostream& os) const override
	{
		os << "GaussSeidel";
	}
};

class ReverseGaussSeidel : public ReverseSOR
{
public: 
	ReverseGaussSeidel() : ReverseSOR(1) {}

	static string Code() { return "rgs"; };

	virtual void Serialize(ostream& os) const override
	{
		os << "Reverse GaussSeidel";
	}
};

class BlockGaussSeidel : public BlockSOR
{
public: 
	BlockGaussSeidel(int blockSize) : BlockSOR(blockSize, 1, Direction::Forward) {}

	static string Code() { return "bgs"; };

	virtual void Serialize(ostream& os) const override
	{
		os << "block-GaussSeidel (blockSize=" << _blockSize << ")";
	}
};

class ReverseBlockGaussSeidel : public BlockSOR
{
public: 
	ReverseBlockGaussSeidel(int blockSize) : BlockSOR(blockSize, 1, Direction::Backward) {}

	static string Code() { return "rbgs"; };

	virtual void Serialize(ostream& os) const override
	{
		os << "Reverse block-GaussSeidel (blockSize=" << _blockSize << ")";
	}
};