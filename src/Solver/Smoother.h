#pragma once
#include <Eigen/Sparse>
#include "IterativeSolver.h"
using namespace std;

class Smoother
{
private:
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

	void Setup(const Eigen::SparseMatrix<double>& A)
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

};