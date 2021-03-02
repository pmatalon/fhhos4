#pragma once
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
		Solver::Setup(A);
		_solver.isSymmetric(true);
		_solver.compute(A);
		this->SetupComputationalWork = Cost::LUFactorization(A);
		Eigen::ComputationInfo info = _solver.info();
		if (info != Eigen::ComputationInfo::Success)
		{
			//cout << "----------------- A -------------------" << A << endl;
			cout << "Error: SparseLU failed to execute with the code " << info << ": " << _solver.lastErrorMessage() << endl;
			exit(EXIT_FAILURE);
		}
	}

	Vector Solve(const Vector& b) override
	{
		Vector x = _solver.solve(b);
		this->SolvingComputationalWork = Cost::LUSolve(*this->Matrix);
		return x;
	}
};