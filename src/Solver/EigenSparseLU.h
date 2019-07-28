#pragma once
#include <Eigen/Sparse>
#include "Solver.h"
using namespace std;

class EigenSparseLU : public Solver
{
private:
	Eigen::SparseLU<SparseMatrix> _solver;

public:
	EigenSparseLU() : Solver() {}

	void Serialize(ostream& os) const override
	{
		os << "Sparse LU factorization (Eigen library)";
	}

	void Setup(const SparseMatrix& A) override
	{
		_solver.compute(A);
		Eigen::ComputationInfo info = _solver.info();
		if (info != Eigen::ComputationInfo::Success)
		{
			cout << "----------------- A -------------------" << A << endl;
			cout << "Error: SparseLU failed to execute with the code " << info << ": " << _solver.lastErrorMessage() << endl;
			exit(EXIT_FAILURE);
		}
	}

	Eigen::VectorXd Solve(const Eigen::VectorXd& b) override
	{
		Eigen::VectorXd x = _solver.solve(b);
		return x;
	}
};