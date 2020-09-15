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
		_solver->StoppingCrit = StoppingCriteria::MaxIterations;
		_solver->MaxIterations = nSmoothingIterations;
		_solver->PrintIterationResults = false;
		_solver->ComputeExactSolution = false;
	}

	int Iterations()
	{
		return _solver->MaxIterations;
	}

	void Setup(const SparseMatrix& A)
	{
		if (_solver->MaxIterations > 0)
			_solver->Setup(A);
	}

	Vector Smooth(Vector& x, const Vector& b)
	{
		return _solver->Solve(b, x);
	}

	friend ostream& operator<<(ostream& os, const Smoother& s)
	{
		os << *(s._solver);
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

class BlockJacobiSmoother : public Smoother
{
public:
	BlockJacobiSmoother(int blockSize, double omega, int nSmoothingIterations) : Smoother(new BlockJacobi(blockSize, omega), nSmoothingIterations) {}
};

class BlockSORSmoother : public Smoother
{
public:
	BlockSORSmoother(int blockSize, double omega, Direction direction, int nSmoothingIterations) : Smoother(new BlockSOR(blockSize, omega, direction), nSmoothingIterations) {}
};