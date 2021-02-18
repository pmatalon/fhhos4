#pragma once
#include "IterativeSolver.h"
#include "Preconditioner.h"
using namespace std;

class ConjugateGradient : public IterativeSolver
{
public:
	Preconditioner Precond;

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
	void Solve(const Vector& b, bool zeroInitialGuess, Vector& initialGuess) override
	{
		const SparseMatrix& A = *this->Matrix;

		this->SolvingComputationalWork = 0; // The work of one iteration of CG is 1 matrix-vector product + 2 dot products

		if (this->ComputeExactSolution)
			this->ExactSolution = this->_directSolver.solve(b);

		IterationResult result = CreateFirstIterationResult(b, initialGuess);

		Vector& x = initialGuess;
		Vector r;
		if (zeroInitialGuess)
			r = b;
		else
		{
			r = b - A * x;                                  result.AddCost(Cost::DAXPY(A));
		}
		result.SetResidual(r);

		double beta = 0;
		Vector z = Precond.Apply(r);                        result.AddCost(Precond.SolvingComputationalWork());
		Vector d = z;
		double r_dot_z = r.dot(z);                          result.AddCost(Cost::Dot(r));
		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			if (this->IterationCount > 0)
				d = z + beta * d;

			Vector Ad = A * d;                                  result.AddCost(Cost::MatVec(A));
			double alpha = r_dot_z / (d.dot(Ad));               result.AddCost(Cost::Dot(r));
			x += alpha * d;

			double old_r_dot_old_z = r_dot_z; // save the dot product before overwriting r and z

			r -= alpha * Ad;
			z = Precond.Apply(r);                               result.AddCost(Precond.SolvingComputationalWork());

			r_dot_z = r.dot(z);                                 result.AddCost(Cost::Dot(r));
			beta = r_dot_z / old_r_dot_old_z;

			this->IterationCount++;

			result.SetX(x);
			result.SetResidual(r);

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();
	}
};