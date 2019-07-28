#pragma once
#include <Eigen/Sparse>
#include "Solver.h"
using namespace std;

class EigenCG : public Solver
{
private:
	Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double>> _solver;

public:
	double Tolerance;

	EigenCG(double tolerance) : Solver() 
	{
		this->Tolerance = tolerance;
	}

	void Serialize(ostream& os) const override
	{
		os << "Conjugate Gradient (Eigen library)";
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

	Eigen::VectorXd Solve(const Eigen::VectorXd& b) override
	{
		Eigen::VectorXd x = _solver.solve(b);
		return x;
	}
};