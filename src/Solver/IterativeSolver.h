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
public:
	Vector ExactSolution;
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

	virtual void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F)
	{
		Utils::FatalError("This solver is not compatible with matrix blocks arising from a hybrid discretization.");
	}

	Vector Solve(const Vector& b) override
	{
		return Solve(b, "0");
	}

	virtual Vector Solve(const Vector& b, string initialGuessCode)
	{
		Vector x;
		bool zeroInitialGuess = false;
		if (initialGuessCode.compare("0") == 0)
		{
			x = Vector::Zero(b.rows());
			zeroInitialGuess = true;
		}
		else if (initialGuessCode.compare("1") == 0)
			x = Vector::Ones(b.rows());
		else if (initialGuessCode.compare("rand") == 0)
			x = Vector::Random(b.rows());
		else
			Utils::FatalError("Unknown initial guess code");
		Solve(b, x, zeroInitialGuess);
		return x;
	}

	virtual void Solve(const Vector& b, Vector& initialGuess, bool zeroInitialGuess)
	{
		const SparseMatrix& A = *this->Matrix;

		this->SolvingComputationalWork = 0;

		if (this->ComputeExactSolution)
			this->ExactSolution = this->_directSolver.solve(b);

		this->IterationCount = 0;
		Vector& x = initialGuess;
		bool xEquals0 = zeroInitialGuess;

		IterationResult result = CreateFirstIterationResult(b, x);
		if (MaxIterations == 0)
			return;

		if (StoppingCrit == StoppingCriteria::NormalizedResidual)
		{
			if (zeroInitialGuess)
				result.SetResidualAsB();
			else
			{
				result.SetResidual(b - A * x);                                  result.AddCost(Cost::DAXPY(A));
			}
		}
		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = ExecuteOneIteration(b, x, xEquals0, result);
			this->IterationCount++;

			if (!result.IsResidualSet() && StoppingCrit == StoppingCriteria::NormalizedResidual)
			{
				result.SetResidual(b - A * x);
				result.AddCost(2 * A.nonZeros());
			}
			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();
	}

	virtual Vector SolveAndComputeResidual(const Vector& b, Vector& initialGuess, bool zeroInitialGuess)
	{
		const SparseMatrix& A = *this->Matrix;

		this->SolvingComputationalWork = 0;

		if (this->ComputeExactSolution)
			this->ExactSolution = this->_directSolver.solve(b);

		this->IterationCount = 0;
		Vector& x = initialGuess;
		bool xEquals0 = zeroInitialGuess;

		Vector r;

		IterationResult result = CreateFirstIterationResult(b, x);
		if (MaxIterations == 0)
		{
			if (zeroInitialGuess)
				r = b;
			else
			{
				r = b - A * x;                                     result.AddCost(Cost::DAXPY(A));
			}
			return r;
		}

		if (MaxIterations == 0)
			return r;

		if (StoppingCrit == StoppingCriteria::NormalizedResidual)
		{
			if (zeroInitialGuess)
			{
				r = b;
				result.SetResidualAsB();
			}
			else
			{
				r = b - A * x;                                     result.AddCost(Cost::DAXPY(A));
				result.SetResidual(r);
			}
		}

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			if (StoppingCrit == StoppingCriteria::MaxIterations && this->IterationCount < this->MaxIterations - 1)
				result = ExecuteOneIteration(b, x, xEquals0, result); // No need to compute the residual until the last iteration
			else
			{
				std::tie(result, r) = ExecuteOneIterationAndComputeResidual(b, x, xEquals0, result);
				if (StoppingCrit == StoppingCriteria::NormalizedResidual)
					result.SetResidual(r);
			}
			this->IterationCount++;

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();
	}

	virtual ~IterativeSolver() {}

protected:
	IterationResult CreateFirstIterationResult(const Vector& b, const Vector& x)
	{
		IterationResult result;
		if (this->MaxIterations == 0 || this->StoppingCrit == StoppingCriteria::MaxIterations)
			return result;
		result.SetA(*this->Matrix);
		result.SetB(b);
		if (this->ComputeExactSolution)
			result.SetExactSolution(this->ExactSolution);
		result.SetX(x);
		result.SetTolerance(this->Tolerance);
		return result;
	}

	virtual IterationResult ExecuteOneIteration(const Vector& b, Vector& x, bool& xEquals0, const IterationResult& oldResult) { assert(false); };
	
	virtual pair<IterationResult, Vector> ExecuteOneIterationAndComputeResidual(const Vector& b, Vector& x, bool& xEquals0, const IterationResult& oldResult)
	{ 
		pair<IterationResult, Vector> p;
		auto&[result, r] = p;

		const SparseMatrix& A = *this->Matrix;

		result = ExecuteOneIteration(b, x, xEquals0, oldResult);
		r = b - A * x;                                         result.AddCost(Cost::DAXPY(A));

		return p;
	}

	bool StoppingCriteriaReached(const IterationResult& result)
	{
		if (MaxIterations == 0)
			return true;
		if (IterationCount == 0)
			return false;
		if (IterationCount >= MaxIterations)
			return true;
		if (StoppingCrit == StoppingCriteria::NormalizedResidual)
		{
			if (std::isinf(result.NormalizedResidualNorm)) // do not remove the prefix std:: or it can become ambiguous according to the compiler
				Utils::FatalError("divergence of the solver");
			if (std::isnan(result.NormalizedResidualNorm))
				Utils::FatalError("the residual is NaN.");
			if (result.NormalizedResidualNorm < this->Tolerance)
				return true;
		}
		return false;
	}
};