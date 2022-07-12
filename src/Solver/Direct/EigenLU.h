#pragma once
#include "../Solver.h"
#include <Eigen/LU>

using namespace std;

class EigenLU : public Solver
{
private:
	Eigen::FullPivLU<DenseMatrix> _solver;

public:
	EigenLU() : Solver() {}

	void Serialize(ostream& os) const override
	{
		os << "FullPivLU factorization (Eigen library)";
	}

	void Setup(const SparseMatrix& A) override
	{
		Utils::FatalError("Solver EigenLU does not take sparse matrices. Use EigenSparseLU instead.");
	}

	void Setup(const DenseMatrix& A) override
	{
		Solver::Setup(A);
		_solver.compute(A);
		//this->SetupComputationalWork = Cost::LUFactorization(A)*1e-6;

		if (_solver.isInvertible())
			Utils::FatalError("Error Eigen::FullPivLU: the matrix is not invertible");
	}

	Vector Solve(const Vector& b) override
	{
		Vector x = _solver.solve(b);
		this->SolvingComputationalWork = 0;//Cost::LUSolve(_solver.matrixL()., _solver.matrixU()); TODO
		return x;
	}
};