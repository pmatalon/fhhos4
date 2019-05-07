#pragma once
#include <Eigen/Sparse>
#include "IterationResult.h"
using namespace std;

class IterativeSolver
{
private:
	Eigen::SparseLU<Eigen::SparseMatrix<double>> _directSolver;
	Eigen::VectorXd _exactSolution;
public:
	double Tolerance = 1e-3;
	int MaxIterations = 200;
	int IterationCount = 0;
	bool PrintIterationResults = true;
	bool ComputeExactSolution = true;

	Eigen::SparseMatrix<double> A;

	IterativeSolver()
	{
	}

	virtual void Setup(const Eigen::SparseMatrix<double>& A)
	{
		this->A = A;
		if (this->ComputeExactSolution)
		{
			this->_directSolver.analyzePattern(this->A);
			this->_directSolver.factorize(this->A);
		}
	}

	Eigen::VectorXd Solve(const Eigen::VectorXd& b)
	{
		Eigen::VectorXd zero = Eigen::VectorXd::Zero(b.rows());
		return Solve(b, zero);
	}

	virtual Eigen::VectorXd Solve(const Eigen::VectorXd& b, Eigen::VectorXd& initialGuess)
	{
		if (this->ComputeExactSolution)
		{
			this->_exactSolution = this->_directSolver.solve(b);
		}

		Eigen::VectorXd x = initialGuess;
		this->IterationCount = 0;

		IterationResult result = SaveIterationResult(x, b);
		while (!StoppingCriteriaReached(result))
		{
			x = ExecuteOneIteration(b, x);
			this->IterationCount++;

			result = SaveIterationResult(x, b);
		}

		return x;
	}

protected:
	IterationResult SaveIterationResult(Eigen::VectorXd& x, const Eigen::VectorXd& b)
	{
		IterationResult result(this->IterationCount, x, b - A * x, b);
		if (this->ComputeExactSolution)
			result.ComputeError(this->_exactSolution);
		if (this->PrintIterationResults)
			cout << result << endl;
		return result;
	}

	virtual Eigen::VectorXd ExecuteOneIteration(const Eigen::VectorXd& b, Eigen::VectorXd& x) = 0;

	bool StoppingCriteriaReached(const IterationResult& result)
	{
		if (IterationCount == 0)
			return false;
		if (IterationCount >= MaxIterations)
			return true;
		if (result.NormalizedResidualNorm < this->Tolerance)
			return true;
		return false;
	}
};