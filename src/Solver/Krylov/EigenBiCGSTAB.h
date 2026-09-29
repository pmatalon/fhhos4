#pragma once
#include "../IterativeSolver.h"
#include <Eigen/IterativeLinearSolvers>
using namespace std;

class EigenBiCGSTAB : public IterativeSolver
{
private:
	Eigen::BiCGSTAB<SparseMatrix, Eigen::DiagonalPreconditioner<double>> _solver;

public:

	EigenBiCGSTAB() : IterativeSolver() {}

	void Serialize(ostream& os) const override
	{
		os << "BiCGSTAB (Eigen library)";
	}

	void Setup(const SparseMatrix& A) override
	{
		Solver::Setup(A);
		_solver.compute(A);
		Eigen::ComputationInfo info = _solver.info();
		if (info != Eigen::ComputationInfo::Success)
		{
			cout << "----------------- A -------------------" << A << endl;
			cout << "Error: ConjugateGradient failed to execute with the code " << info << "." << endl;
			exit(EXIT_FAILURE);
		}
		_solver.setTolerance(this->Tolerance);
		_solver.setMaxIterations(this->MaxIterations);
	}

	void Solve(const Vector& b, Vector& x, bool xEquals0, bool computeResidual, bool computeAx) override
	{
#if EIGEN_VERSION_AT_LEAST(3, 5, 0) // Eigen 5 reports itself as 3.5 in these macros
		// Eigen 5's BiCGSTAB compares the tolerance to the absolute residual norm (Eigen 3.x used
		// ||r|| <= tol * ||b||). Scale it to keep a relative criterion: otherwise, a right-hand side
		// already smaller than tol makes it return x = 0 immediately (e.g. when used as an inner
		// preconditioner solver near convergence, which then breaks the outer Krylov method).
		_solver.setTolerance(this->Tolerance * b.norm());
#endif
		if (xEquals0)
			x = _solver.solve(b);
		else
			x = _solver.solveWithGuess(b, x);
		this->IterationCount = _solver.iterations();
	}

	Vector Solve(const Vector& b) override
	{
		return IterativeSolver::Solve(b);
	}
};