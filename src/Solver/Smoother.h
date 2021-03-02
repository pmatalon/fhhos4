#pragma once
#include "BlockSOR.h"
#include "BlockGaussSeidel.h"
#include "BlockJacobi.h"
#include "GaussSeidel.h"
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

	bool CanOptimizeResidualComputation()
	{
		return _solver->CanOptimizeResidualComputation() && _solver->MaxIterations > 0;
	}

	void Smooth(Vector& x, const Vector& b, bool& xEquals0)
	{
		_solver->Solve(b, x, xEquals0, false, false);
		xEquals0 = false;
	}

	void SmoothAndComputeResidualOrAx(Vector& x, const Vector& b, bool& xEquals0, bool computeResidual, bool computeAx)
	{
		_solver->Solve(b, x, xEquals0, computeResidual, computeAx);
		xEquals0 = false;
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

	Vector&& Residual()
	{
		return std::move(_solver->Residual);
	}

	Vector&& Ax()
	{
		return std::move(_solver->Ax);
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

class BlockGaussSeidelSmoother : public Smoother
{
public:
	BlockGaussSeidelSmoother(int blockSize, Direction direction, int nSmoothingIterations) : Smoother(new BlockGaussSeidel(blockSize, direction), nSmoothingIterations) {}
};

class GaussSeidelSmoother : public Smoother
{
public:
	GaussSeidelSmoother(Direction direction, int nSmoothingIterations) : Smoother(new GaussSeidel(direction), nSmoothingIterations) {}
};