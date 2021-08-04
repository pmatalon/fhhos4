#pragma once
#include "../Solver.h"
#include <Eigen/SparseCholesky>
using namespace std;

class EigenSparseCholesky : public Solver
{
private:
	Eigen::SimplicialLDLT<ColMajorSparseMatrix> _solver; // Eigen only supports ColMajor so far...

public:
	EigenSparseCholesky() : Solver() {}

	void Serialize(ostream& os) const override
	{
		os << "Cholesky factorization (Eigen library)";
	}

	void Setup(const SparseMatrix& A) override
	{
		Solver::Setup(A);
		if (A.IsRowMajor)
			_solver.compute(A.transpose());
		else
			_solver.compute(A);
		this->SetupComputationalWork = Cost::CholeskyFactorization(A);
		Eigen::ComputationInfo info = _solver.info();
		if (info != Eigen::ComputationInfo::Success)
		{
			//cout << "----------------- A -------------------" << A << endl;
			cout << "Error: SimplicialLDLT failed to execute with the code " << info << endl;
			exit(EXIT_FAILURE);
		}
	}

	Vector Solve(const Vector& b) override
	{
		Vector x = _solver.solve(b);
		this->SolvingComputationalWork = Cost::CholeskySolve(_solver.matrixL());
		return x;
	}
};