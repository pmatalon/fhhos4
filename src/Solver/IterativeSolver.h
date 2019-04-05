#pragma once
#include <Eigen/Sparse>
using namespace std;

class IterativeSolver
{
private:
public:
	double Tolerance = 1e-3;
	int MaxIterations = 200;
	int IterationCount = 0;

	IterativeSolver()
	{
	}

	virtual void Setup(const Eigen::SparseMatrix<double>& A) = 0;

	virtual Eigen::VectorXd Solve(const Eigen::VectorXd& b, Eigen::VectorXd& initialGuess)
	{
		Eigen::VectorXd x = initialGuess;
		this->IterationCount = 0;

		while (!StoppingCriteriaReached())
		{
			x = ExecuteOneIteration(b, x);
			this->IterationCount++;
		}

		return x;
	}

protected:
	virtual Eigen::VectorXd ExecuteOneIteration(const Eigen::VectorXd& b, Eigen::VectorXd& x) = 0;

	bool StoppingCriteriaReached()
	{
		if (IterationCount > MaxIterations)
			return true;
	}
};