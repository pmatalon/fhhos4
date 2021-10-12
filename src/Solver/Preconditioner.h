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

	bool CanOptimizeResidualComputation()
	{
		return !_solver ? false : _solver->CanOptimizeResidualComputation();
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

	Vector Apply(const Vector& r)
	{
		if (_solver)
			return _solver->Solve(r);
		return r;
	}

	pair<Vector, Vector> ApplyAndComputeAe(const Vector& r)
	{
		assert(_solver);

		pair<Vector, Vector> p;
		auto&[e, Ae] = p;

		e = Vector::Zero(r.rows());
		bool eEquals0 = true;
		_solver->Solve(r, e, eEquals0, false, true);
		Ae = std::move(_solver->Ax);

		return p;
	}
	
	MFlops SetupComputationalWork()
	{
		if (_solver)
			return _solver->SetupComputationalWork;
		return 0;
	}

	MFlops SolvingComputationalWork()
	{
		if (_solver)
			return _solver->SolvingComputationalWork;
		return 0;
	}
};