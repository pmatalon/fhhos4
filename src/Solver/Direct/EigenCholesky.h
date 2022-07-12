#pragma once
#include "../Solver.h"

using namespace std;

class EigenCholesky : public Solver
{
private:
	Eigen::LLT<DenseMatrix> _solver;

public:
	EigenCholesky() : Solver() {}

	void Serialize(ostream& os) const override
	{
		os << "Cholesky factorization (Eigen library)";
	}

	void Setup(const SparseMatrix& A) override
	{
		Utils::FatalError("Solver EigenCholesky does not take sparse matrices. Use EigenSparseCholesky instead.");
	}

	void Setup(const DenseMatrix& A) override
	{
		Solver::Setup(A);
		_solver.compute(A);
		//this->SetupComputationalWork = Cost::LUFactorization(A)*1e-6;

		Eigen::ComputationInfo info = _solver.info();
		if (info != Eigen::ComputationInfo::Success)
		{
			//cout << "----------------- A -------------------" << A << endl;
			cout << "Error: Eigen::LLT failed to execute with the code " << info << endl;
			exit(EXIT_FAILURE);
		}
	}

	Vector Solve(const Vector& b) override
	{
		Vector x = _solver.solve(b);
		this->SolvingComputationalWork = 0;//Cost::LUSolve(_solver.matrixL()., _solver.matrixU()); TODO
		return x;
	}
};