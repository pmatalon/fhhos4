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
		_solver.compute(A);
		this->SetupComputationalWork = 2.0 / 3.0 * pow(A.rows(), 3); // LU factorization
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
		this->SolvingComputationalWork = 0;
		Vector x = _solver.solve(b);
		this->SolvingComputationalWork = b.rows() * b.rows(); // back substitution
		return x;
	}
};