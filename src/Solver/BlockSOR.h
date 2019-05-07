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
		IterativeSolver::Setup(A);

		D = Eigen::SparseMatrix<double>(A.rows(), A.cols());
		L = Eigen::SparseMatrix<double>(A.rows(), A.cols());
		U = Eigen::SparseMatrix<double>(A.rows(), A.cols());

		NonZeroCoefficients D_coeffs(A.nonZeros());
		NonZeroCoefficients L_coeffs(A.nonZeros());
		NonZeroCoefficients U_coeffs(A.nonZeros());

		for (unsigned int i = 0; i < A.rows() / _blockSize; ++i)
		{
			for (unsigned int j = 0; j < A.rows() / _blockSize; ++j)
			{
				unsigned int k = i * _blockSize;
				unsigned int l = j * _blockSize;
				auto block = A.block(k, l, _blockSize, _blockSize);
				if (i == j)
					D_coeffs.Add(k, l, block);
				else if (i < j)
					U_coeffs.Add(k, l, block);
				else
					L_coeffs.Add(k, l, block);
			}
		}

		D_coeffs.Fill(D);
		U_coeffs.Fill(U);
		L_coeffs.Fill(L);

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

class BlockGaussSeidel : public BlockSOR
{
public: BlockGaussSeidel(int blockSize) : BlockSOR(blockSize, 1, Direction::Forward) {}
};

class ReverseBlockGaussSeidel : public BlockSOR
{
public: ReverseBlockGaussSeidel(int blockSize) : BlockSOR(blockSize, 1, Direction::Backward) {}
};