#pragma once
#include "IterativeSolver.h"
using namespace std;

class Preconditioner
{
protected:
	IterativeSolver* _solver = nullptr;
public:
	Preconditioner() {}

	Preconditioner(IterativeSolver* solver)
	{
		_solver = solver;
		_solver->ComputeExactSolution = false;
		_solver->PrintIterationResults = false;
		_solver->StoppingCrit = StoppingCriteria::MaxIterations;
		_solver->MaxIterations = 1;
	}

	IterativeSolver* GetSolver()
	{
		return _solver;
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

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F)
	{
		if (_solver)
			_solver->Setup(A, A_T_T, A_T_F, A_F_F);
	}

	Vector Solve(const Vector& b)
	{
		if (_solver)
			return _solver->Solve(b);
		return b;
	}
	
	BigNumber SetupComputationalWork()
	{
		if (_solver)
			return _solver->SetupComputationalWork;
		return 0;
	}

	BigNumber SolvingComputationalWork()
	{
		if (_solver)
			return _solver->SolvingComputationalWork;
		return 0;
	}
};