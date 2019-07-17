#pragma once
#include <Eigen/Sparse>
using namespace std;

struct IterationResult
{
	int IterationNumber;
	Eigen::VectorXd x;
	Eigen::VectorXd r;
	Eigen::VectorXd e;
	double ResidualNorm;
	double NormalizedResidualNorm;
	double RelativeErrorNorm = -1;

	IterationResult(int iterationNumber, const Eigen::VectorXd& x, const Eigen::VectorXd& r, const Eigen::VectorXd& b)
	{
		this->IterationNumber = iterationNumber;
		this->x = x;
		this->r = r;
		this->ResidualNorm = r.norm();
		double b_norm = b.norm();
		this->NormalizedResidualNorm = b_norm > 0 ? ResidualNorm / b.norm() : ResidualNorm;
	}

	void ComputeError(const Eigen::VectorXd& exactSolution)
	{
		this->e = exactSolution - x;
		double exactSolNorm = exactSolution.norm();
		this->RelativeErrorNorm = exactSolNorm > 0 ? e.norm() / exactSolution.norm() : e.norm();
	}

	friend ostream& operator<<(ostream& os, const IterationResult& result)
	{
		if (result.IterationNumber == 0)
		{
			os << "It.\tNormalized res";
			if (result.RelativeErrorNorm != -1)
				os << "\tRelative err";
		}
		else
		{
			os << result.IterationNumber << "\t" << std::scientific << result.NormalizedResidualNorm;
			if (result.RelativeErrorNorm != -1)
				os << "\t" << result.RelativeErrorNorm;
		}
		return os;
	}
};