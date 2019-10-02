#pragma once
#include "BlockSOR.h"
#include "BlockJacobi.h"
using namespace std;

class Smoother
{
protected:
	IterativeSolver* _solver;
	int _nSmoothingIterations;

public:

	Smoother(IterativeSolver* solver, int nSmoothingIterations) :
		_solver(solver)
	{
		_solver->MaxIterations = nSmoothingIterations;
		_solver->PrintIterationResults = false;
		_solver->ComputeExactSolution = false;
	}

	void Setup(const SparseMatrix& A)
	{
		_solver->Setup(A);
	}

	Eigen::VectorXd Smooth(Eigen::VectorXd& x, const Eigen::VectorXd& b)
	{
		return _solver->Solve(b, x);
	}

	friend ostream& operator<<(ostream& os, const Smoother& s)
	{
		os << s._solver->MaxIterations << " sweep(s) of " << *(s._solver);
		return os;
	}

	BigNumber SolvingComputationalWork()
	{
		return _solver->SolvingComputationalWork;
	}

	virtual ~Smoother()
	{
		delete _solver;
	}

};

class JacobiSmoother : public Smoother
{
public:
	JacobiSmoother(int nSmoothingIterations) : Smoother(new Jacobi(), nSmoothingIterations) {}
};

class BlockJacobiSmoother : public Smoother
{
public:
	BlockJacobiSmoother(int blockSize, int nSmoothingIterations) : Smoother(new BlockJacobi(blockSize), nSmoothingIterations) {}
};

class GaussSeidelSmoother : public Smoother
{
public:
	GaussSeidelSmoother(int nSmoothingIterations) : Smoother(new GaussSeidel(), nSmoothingIterations) {}
};

class ReverseGaussSeidelSmoother : public Smoother
{
public:
	ReverseGaussSeidelSmoother(int nSmoothingIterations) : Smoother(new ReverseGaussSeidel(), nSmoothingIterations) {}
};

class BlockGaussSeidelSmoother : public Smoother
{
public:
	BlockGaussSeidelSmoother(int blockSize, int nSmoothingIterations) : Smoother(new BlockGaussSeidel(blockSize), nSmoothingIterations) {}
};

class ReverseBlockGaussSeidelSmoother : public Smoother
{
public:
	ReverseBlockGaussSeidelSmoother(int blockSize, int nSmoothingIterations) : Smoother(new ReverseBlockGaussSeidel(blockSize), nSmoothingIterations) {}
};
