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
protected:
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

	virtual void Serialize(ostream& os) const override
	{
		os << "Block-SOR (blockSize=" << _blockSize << ", omega=" << _omega << ", direction=" << (_direction == Direction::Forward ? "forward" : "backward") << ")";
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

		for (int k = 0; k < A.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
			{
				auto iBlock = it.row() / _blockSize;
				auto jBlock = it.col() / _blockSize;
				if (iBlock == jBlock)
					D_coeffs.Add(it.row(), it.col(), it.value());
				else if (iBlock < jBlock)
					U_coeffs.Add(it.row(), it.col(), it.value());
				else
					L_coeffs.Add(it.row(), it.col(), it.value());
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

	virtual void Serialize(ostream& os) const override
	{
		os << "GaussSeidel";
	}
};

class ReverseGaussSeidel : public ReverseSOR
{
public: 
	ReverseGaussSeidel() : ReverseSOR(1) {}

	virtual void Serialize(ostream& os) const override
	{
		os << "Reverse GaussSeidel";
	}
};

class BlockGaussSeidel : public BlockSOR
{
public: 
	BlockGaussSeidel(int blockSize) : BlockSOR(blockSize, 1, Direction::Forward) {}

	virtual void Serialize(ostream& os) const override
	{
		os << "block-GaussSeidel (blockSize=" << _blockSize << ")";
	}
};

class ReverseBlockGaussSeidel : public BlockSOR
{
public: 
	ReverseBlockGaussSeidel(int blockSize) : BlockSOR(blockSize, 1, Direction::Backward) {}

	virtual void Serialize(ostream& os) const override
	{
		os << "Reverse block-GaussSeidel (blockSize=" << _blockSize << ")";
	}
};