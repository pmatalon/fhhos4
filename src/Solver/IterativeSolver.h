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
	bool ComputeExactSolution = false;

	SparseMatrix A;

	IterativeSolver() : Solver() {}

	virtual void Setup(const SparseMatrix& A) override
	{
		this->A = A;
		if (this->ComputeExactSolution)
			this->_directSolver.compute(A);
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
		this->SolvingComputationalWork = 0;

		if (this->ComputeExactSolution)
			this->_exactSolution = this->_directSolver.solve(b);

		this->IterationCount = 0;

		IterationResult result = CreateFirstIterationResult(b, initialGuess);
		result.SetResidual(b - A * initialGuess);
		result.AddCost(2 * A.nonZeros());
		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			auto x = result.X();
			result = ExecuteOneIteration(b, x, result);
			this->IterationCount++;

			if (!result.IsResidualSet())
			{
				result.SetResidual(b - A * result.X());
				result.AddCost(2 * A.nonZeros());
			}
			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();

		return result.X();
	}

	virtual ~IterativeSolver() {}

protected:
	IterationResult CreateFirstIterationResult(const Eigen::VectorXd& b, const Eigen::VectorXd& x)
	{
		IterationResult result;
		result.SetB(b);
		if (this->ComputeExactSolution)
			result.SetExactSolution(this->_exactSolution);
		result.SetX(x);
		return result;
	}
	/*IterationResult SaveIterationResult(Eigen::VectorXd& x, const Eigen::VectorXd& b, const Eigen::VectorXd& r)
	{
		IterationResult result = CreateNewIterationResult(x, b, r);
		if (this->PrintIterationResults)
			cout << result << endl;
		return result;
	}

	IterationResult SaveIterationResult(Eigen::VectorXd& x, const Eigen::VectorXd& b)
	{
		return SaveIterationResult(x, b, b - A * x);
	}*/

	virtual IterationResult ExecuteOneIteration(const Eigen::VectorXd& b, Eigen::VectorXd& x, const IterationResult& oldResult) {};

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