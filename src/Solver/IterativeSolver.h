#pragma once
#include <Eigen/Sparse>
#include "Solver.h"
#include "IterationResult.h"
using namespace std;

class IterativeSolver : public Solver
{
protected:
	Eigen::SparseLU<SparseMatrix> _directSolver;
	Eigen::VectorXd _exactSolution;
public:
	double Tolerance = 1e-3;
	int MaxIterations = 200;
	int IterationCount = 0;
	bool PrintIterationResults = true;
	bool ComputeExactSolution = true;

	SparseMatrix A;

	IterativeSolver() : Solver() {}

	virtual void Setup(const SparseMatrix& A) override
	{
		this->A = A;
	}

	Eigen::VectorXd Solve(const Eigen::VectorXd& b) override
	{
		return Solve(b, "0");
	}

	Eigen::VectorXd Solve(const Eigen::VectorXd& b, string initialGuessCode)
	{
		Eigen::VectorXd initialGuess;
		if (initialGuessCode.compare("0") == 0)
			initialGuess = Eigen::VectorXd::Zero(b.rows());
		else if (initialGuessCode.compare("1") == 0)
			initialGuess = Eigen::VectorXd::Ones(b.rows());
		else
			assert(false);
		return Solve(b, initialGuess);
	}

	virtual Eigen::VectorXd Solve(const Eigen::VectorXd& b, Eigen::VectorXd& initialGuess)
	{
		if (this->ComputeExactSolution)
			this->_exactSolution = this->_directSolver.solve(b);

		Eigen::VectorXd x = initialGuess;
		this->IterationCount = 0;

		IterationResult result = SaveIterationResult(x, b);
		while (!StoppingCriteriaReached(result))
		{
			x = ExecuteOneIteration(b, x);
			this->IterationCount++;

			result = SaveIterationResult(x, b);
		}

		if (this->PrintIterationResults)
			cout << endl;

		return x;
	}

	virtual ~IterativeSolver() {}

protected:
	IterationResult SaveIterationResult(Eigen::VectorXd& x, const Eigen::VectorXd& b, const Eigen::VectorXd& r)
	{
		IterationResult result(this->IterationCount, x, r, b);
		if (this->ComputeExactSolution)
			result.ComputeError(this->_exactSolution);
		if (this->PrintIterationResults)
			cout << result << endl;
		return result;
	}

	IterationResult SaveIterationResult(Eigen::VectorXd& x, const Eigen::VectorXd& b)
	{
		return SaveIterationResult(x, b, b - A * x);
	}

	virtual Eigen::VectorXd ExecuteOneIteration(const Eigen::VectorXd& b, Eigen::VectorXd& x) {};

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