#pragma once
#include "../IterativeSolver.h"
#include "../Preconditioner.h"
using namespace std;

class ConjugateGradient : public IterativeSolver
{
public:
	SolverPreconditioner Precond;

	ConjugateGradient()
	{ }

	virtual void Serialize(ostream& os) const override
	{
		os << "Conjugate Gradient, preconditioner: " << Precond;
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		this->Precond.Setup(A);
		this->SetupComputationalWork = this->Precond.SetupComputationalWork();
	}

private:
	void Solve(const Vector& b, Vector& initialGuess, bool zeroInitialGuess) override
	{
		const SparseMatrix& A = *this->Matrix;

		this->SolvingComputationalWork = 0; // The work of one iteration of CG is 1 matrix-vector product + 2 dot products

		if (this->ComputeExactSolution)
			this->ExactSolution = this->_directSolver.solve(b);

		IterationResult result = CreateFirstIterationResult(b, initialGuess);

		Vector& x = initialGuess;
		Vector r;
		if (zeroInitialGuess)
		{
			r = b;
			result.SetResidualAsB();
		}
		else
		{
			r = b - A * x;                                  result.AddWorkInFlops(Cost::DAXPY(A));
			result.SetResidualNorm(r.norm());               result.AddWorkInFlops(Cost::Norm(r));
		}

		double beta = 0;
		Vector z = Precond.Apply(r);                        result.AddWorkInMFlops(Precond.SolvingComputationalWork());
		Vector d = z;
		double r_dot_z = r.dot(z);                          result.AddWorkInFlops(Cost::Dot(r));
		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			if (this->IterationCount > 0)
			{
				// Update direction of research
				d = z + beta * d;                               // TODO AddWork
			}

			Vector Ad = A * d;                                  result.AddWorkInFlops(Cost::MatVec(A));

			// Step length in the direction of research
			double alpha = r_dot_z / (d.dot(Ad));               result.AddWorkInFlops(Cost::Dot(r));

			// Moving from the current solution to the next
			// by taking the step in the direction of research
			x += alpha * d;                                     result.AddWorkInFlops(Cost::VectorDAXPY(x));

			double old_r_dot_old_z = r_dot_z; // save the dot product before overwriting r and z

			// New residual
			r -= alpha * Ad;                                    result.AddWorkInFlops(Cost::VectorDAXPY(r));
			z = Precond.Apply(r);                               result.AddWorkInMFlops(Precond.SolvingComputationalWork());

			// New dot product
			r_dot_z = r.dot(z);                                 result.AddWorkInFlops(Cost::Dot(r));
			// New step for the research direction
			beta = r_dot_z / old_r_dot_old_z;

			this->IterationCount++;

			result.SetX(x);
			result.SetResidualNorm(r.norm());                   result.AddWorkInFlops(Cost::Norm(r));

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();
	}
};