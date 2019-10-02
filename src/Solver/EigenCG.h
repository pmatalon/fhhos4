#pragma once
#include <Eigen/Sparse>
#include "IterativeSolver.h"
using namespace std;

class EigenCG : public IterativeSolver
{
private:
	Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double>> _solver;

public:

	EigenCG() : IterativeSolver() {}

	void Serialize(ostream& os) const override
	{
		os << "Conjugate Gradient (Eigen library), diagonal preconditioner";
	}

	void Setup(const SparseMatrix& A) override
	{
		_solver.compute(A);
		Eigen::ComputationInfo info = _solver.info();
		if (info != Eigen::ComputationInfo::Success)
		{
			cout << "----------------- A -------------------" << A << endl;
			cout << "Error: ConjugateGradient failed to execute with the code " << info << "." << endl;
			exit(EXIT_FAILURE);
		}
		_solver.setTolerance(this->Tolerance);
	}

	Vector Solve(const Vector& b, Vector& initialGuess) override
	{
		Vector x = _solver.solveWithGuess(b, initialGuess);
		this->IterationCount = _solver.iterations();
		return x;
	}
};