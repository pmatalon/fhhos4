#pragma once
#include <Eigen/Sparse>
#include "IterativeSolver.h"
using namespace std;

enum class Direction : unsigned
{
	Forward = 0,
	Backward = 1 << 1
};

class BlockSOR : public IterativeSolver
{
private:
	int _blockSize;
	double _omega;
	Direction _direction;

	Eigen::SparseMatrix<double> D;
	Eigen::SparseMatrix<double> L;
	Eigen::SparseMatrix<double> U;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> _solver;
public:

	BlockSOR(int blockSize, double omega, Direction direction)
	{
		this->_blockSize = blockSize;
		this->_omega = omega;
		this->_direction = direction;
	}

	void Setup(const Eigen::SparseMatrix<double>& A) override
	{
		D = Eigen::SparseMatrix<double>(A.rows(), A.cols());
		L = Eigen::SparseMatrix<double>(A.rows(), A.cols());
		U = Eigen::SparseMatrix<double>(A.rows(), A.cols());

		for (unsigned int i = 0; i < A.rows() / _blockSize; ++i)
		{
			for (unsigned int j = 0; j < A.rows() / _blockSize; ++j)
			{
				unsigned int k = i * _blockSize;
				unsigned int l = j * _blockSize;
				auto block = A.block(k, l, _blockSize, _blockSize);
				/*if (i == j)
					D.block(k, l, _blockSize, _blockSize) = block;
				else if (i < j)
					U.block(k, l, _blockSize, _blockSize) = block;
				else
					L.block(k, l,_ blockSize, _blockSize) = block;*/
			}
		}

		Eigen::SparseMatrix<double> M = _direction == Direction::Forward ? (D + _omega * L) : (D + _omega * U);
		_solver.analyzePattern(M);
		_solver.factorize(M);
	}

private:
	Eigen::VectorXd ExecuteOneIteration(const Eigen::VectorXd& b, Eigen::VectorXd& x) override
	{
		Eigen::VectorXd rhs;
		if (_direction == Direction::Forward)
			rhs = (-_omega * U + (1 - _omega)*D)*x + _omega * b;
		else
			rhs = (-_omega * L + (1 - _omega)*D)*x + _omega * b;

		x = _solver.solve(rhs);
		return x;
	}
};

class SOR : public BlockSOR
{
public: SOR(double omega) : BlockSOR(1, omega, Direction::Forward) {}
};

class ReverseSOR : public BlockSOR
{
public: ReverseSOR(double omega) : BlockSOR(1, omega, Direction::Backward) {}
};

class GaussSeidel : public SOR
{
public: GaussSeidel() : SOR(1) {}
};

class ReverseGaussSeidel : public ReverseSOR
{
public: ReverseGaussSeidel() : ReverseSOR(1) {}
};