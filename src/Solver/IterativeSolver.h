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
	Vector Residual;
	Vector Ax;

	IterativeSolver() : Solver() {}

	virtual void Setup(const SparseMatrix& A) override
	{
		Solver::Setup(A);
		if (this->ComputeExactSolution)
			this->_directSolver.compute(A);
	}

	virtual void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F)
	{
		Utils::FatalError("This solver is not compatible with matrix blocks arising from a hybrid discretization.");
	}

	virtual bool CanOptimizeResidualComputation()
	{
		return false;
	}

	Vector Solve() override
	{
		Utils::FatalError("Solve() not implemented in this solver.");
		return Vector::Zero(0);
	}

	Vector Solve(const Vector& b) override
	{
		return Solve(b, "0");
	}

	virtual Vector Solve(const Vector& b, string initialGuessCode)
	{
		Vector x;
		bool xEquals0 = false;
		if (initialGuessCode.compare("0") == 0)
		{
			x = Vector::Zero(b.rows());
			xEquals0 = true;
		}
		else if (initialGuessCode.compare("1") == 0)
			x = Vector::Ones(b.rows());
		else if (initialGuessCode.compare("rand") == 0)
			x = Vector::Random(b.rows());
		else
			Utils::FatalError("Unknown initial guess code");

		Solve(b, x, xEquals0);
		return x;
	}

	virtual void Solve(const Vector& b, Vector& x, bool xEquals0)
	{
		bool computeResidual = false;
		bool computeAx = false;
		Solve(b, x, xEquals0, computeResidual, computeAx);
	}

	virtual void Solve(const Vector& b, Vector& x, bool xEquals0, bool computeResidual, bool computeAx)
	{
		DoBeforeSolving();

		const SparseMatrix& A = *this->Matrix;

		this->SolvingComputationalWork = 0;
		this->Residual = Vector();

		if (this->ComputeExactSolution)
			this->ExactSolution = this->_directSolver.solve(b);

		this->IterationCount = 0;
		
		IterationResult result = CreateFirstIterationResult(b, x);
		if (MaxIterations == 0)
		{
			if (computeResidual && !computeAx)
			{
				if (xEquals0)
					this->Residual = b;
				else
				{
					this->Residual = b - A * x;                              result.AddWorkInFlops(Cost::DAXPY(A));
				}
			}
			else if (computeResidual && computeAx)
			{
				if (xEquals0)
				{
					this->Residual = b;
					this->Ax = Vector::Zero(b.rows());
				}
				else
				{
					this->Ax = A * x;                                        result.AddWorkInFlops(Cost::MatVec(A));
					this->Residual = b - this->Ax;                           result.AddWorkInFlops(Cost::AddVec(b));
				}
			}
			else if (computeAx)
			{
				if (xEquals0)
					this->Ax = Vector::Zero(b.rows());
				else
				{
					this->Ax = A * x;                                        result.AddWorkInFlops(Cost::MatVec(A));
				}
			}
			return;
		}

		if (StoppingCrit == StoppingCriteria::NormalizedResidual)
		{
			if (xEquals0)
				result.SetResidualAsB();
			else
			{
				this->Residual = b - A * x;                                          result.AddWorkInFlops(Cost::DAXPY(A));
				result.SetResidualNorm(this->Residual.norm());                       result.AddWorkInFlops(Cost::Norm(b));
			}
		}

		if (this->PrintIterationResults)
			cout << result << endl;


		while (!StoppingCriteriaReached(result))
		{
			if (!CanOptimizeResidualComputation())
			{
				result = ExecuteOneIteration(b, x, xEquals0, false, false, result);
				if (!result.IsResidualSet() && StoppingCrit == StoppingCriteria::NormalizedResidual)
				{
					this->Residual = b - A * x;                                          result.AddWorkInFlops(Cost::DAXPY(A));
					result.SetResidualNorm(this->Residual.norm());                       result.AddWorkInFlops(Cost::Norm(b));
				}                                      
			}
			else
			{
				if (StoppingCrit == StoppingCriteria::MaxIterations)
				{
					if ((!computeResidual && !computeAx) || this->IterationCount < this->MaxIterations - 1)
						result = ExecuteOneIteration(b, x, xEquals0, false, false, result); // No need to compute the residual until the last iteration
					else
						result = ExecuteOneIteration(b, x, xEquals0, computeResidual, computeAx, result);
				}
				else
				{
					result = ExecuteOneIteration(b, x, xEquals0, true, computeAx, result);
					assert(result.Residual.rows() > 0);
					result.SetResidualNorm(result.Residual.norm());                                                result.AddWorkInFlops(Cost::Norm(b));
				}
			}
			this->IterationCount++;

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();
		if (computeResidual)
			this->Residual = std::move(result.Residual);
		if (computeAx)
			this->Ax = std::move(result.Ax);
	}

	virtual ~IterativeSolver() {}

protected:
	virtual void DoBeforeSolving()
	{}

	IterationResult CreateFirstIterationResult(const Vector& b, const Vector& x)
	{
		IterationResult result;
		if (this->MaxIterations == 0 || this->StoppingCrit == StoppingCriteria::MaxIterations)
			return result;
		if (this->Matrix)
			result.SetA(*this->Matrix);
		result.SetB(b);
		if (this->ComputeExactSolution)
			result.SetExactSolution(this->ExactSolution);
		result.SetX(x);
		result.SetTolerance(this->Tolerance);
		return result;
	}

	virtual IterationResult ExecuteOneIteration(const Vector& b, Vector& x, bool& xEquals0, bool computeResidual, bool computeAx, const IterationResult& oldResult) 
	{ 
		Utils::FatalError("The function Solve() or ExecuteOneIteration() must be overriden in the subclass.");
		return IterationResult(); // to avoid warning
	};


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