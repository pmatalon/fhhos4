#pragma once
#include "Solver.h"
#include "IterationResult.h"
using namespace std;

enum class StoppingCriteria : unsigned
{
	NormalizedResidual = 0,
	MaxIterations = 1
};

class IterativeSolver : public Solver
{
protected:
	Eigen::SparseLU<SparseMatrix> _directSolver;
	Vector _exactSolution;
public:
	double Tolerance = 1e-3;
	int MaxIterations = 200;
	int IterationCount = 0;
	bool PrintIterationResults = true;
	bool ComputeExactSolution = false;
	StoppingCriteria StoppingCrit = StoppingCriteria::NormalizedResidual;

	const SparseMatrix* Matrix;

	IterativeSolver() : Solver() {}

	virtual void Setup(const SparseMatrix& A) override
	{
		this->Matrix = &A;
		if (this->ComputeExactSolution)
			this->_directSolver.compute(A);
	}

	Vector Solve(const Vector& b) override
	{
		return Solve(b, "0");
	}

	virtual Vector Solve(const Vector& b, string initialGuessCode)
	{
		Vector initialGuess;
		if (initialGuessCode.compare("0") == 0)
			initialGuess = Vector::Zero(b.rows());
		else if (initialGuessCode.compare("1") == 0)
			initialGuess = Vector::Ones(b.rows());
		else if (initialGuessCode.compare("rand") == 0)
			initialGuess = Vector::Random(b.rows());
		else
			assert(false);
		return Solve(b, initialGuess);
	}

	virtual Vector Solve(const Vector& b, Vector& initialGuess)
	{
		const SparseMatrix& A = *this->Matrix;

		this->SolvingComputationalWork = 0;

		if (this->ComputeExactSolution)
			this->_exactSolution = this->_directSolver.solve(b);

		this->IterationCount = 0;

		IterationResult result = CreateFirstIterationResult(b, initialGuess);
		if (MaxIterations == 0)
		{
			result.SetX(initialGuess);
			return initialGuess;
		}

		if (StoppingCrit == StoppingCriteria::NormalizedResidual)
		{
			result.SetResidual(b - A * initialGuess);
			result.AddCost(2 * A.nonZeros());
		}
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
	IterationResult CreateFirstIterationResult(const Vector& b, const Vector& x)
	{
		IterationResult result;
		result.SetB(b);
		if (this->ComputeExactSolution)
			result.SetExactSolution(this->_exactSolution);
		result.SetX(x);
		result.SetTolerance(this->Tolerance);
		return result;
	}

	virtual IterationResult ExecuteOneIteration(const Vector& b, Vector& x, const IterationResult& oldResult) { assert(false); };

	bool StoppingCriteriaReached(const IterationResult& result)
	{
		if (MaxIterations == 0)
			return true;
		if (IterationCount == 0)
			return false;
		if (IterationCount >= MaxIterations)
			return true;
		if (std::isinf(result.NormalizedResidualNorm)) // do not remove the prefix std:: or it can become ambiguous according to the compiler
			Utils::FatalError("divergence of the solver");
		if (StoppingCrit == StoppingCriteria::NormalizedResidual && result.NormalizedResidualNorm < this->Tolerance)
			return true;
		return false;
	}
};