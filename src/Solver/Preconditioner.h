#pragma once
#include "IterativeSolver.h"
using namespace std;

class Preconditioner
{
private:
	IterativeSolver* _solver = nullptr;
public:
	Preconditioner() {}

	Preconditioner(IterativeSolver* solver)
	{
		_solver = solver;
		_solver->ComputeExactSolution = false;
		_solver->PrintIterationResults = false;
		_solver->MaxIterations = 1;
	}

	bool IsIdentity()
	{
		return _solver == nullptr;
	}

	friend ostream& operator<<(ostream& os, const Preconditioner& p)
	{
		if (p._solver)
			os << *(p._solver);
		else
			os << "none";
		return os;
	}

	void Setup(const SparseMatrix& A)
	{
		if (_solver)
			_solver->Setup(A);
	}

	Eigen::VectorXd Solve(const Eigen::VectorXd& b)
	{
		if (_solver)
			return _solver->Solve(b);
		return b;
	}
	
	BigNumber SetupComputationalWork()
	{
		return _solver->SetupComputationalWork;
	}

	BigNumber SolvingComputationalWork()
	{
		return _solver->SolvingComputationalWork;
	}
};