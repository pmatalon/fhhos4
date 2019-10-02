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
	Vector Solve(const Vector& b, Vector& initialGuess) override
	{
		this->SolvingComputationalWork = 0;

		if (this->ComputeExactSolution)
			this->_exactSolution = this->_directSolver.solve(b);

		IterationResult result = CreateFirstIterationResult(b, initialGuess);

		Vector x = initialGuess;
		Vector r = b - A * x;
		result.AddCost(2 * A.nonZeros());
		result.SetResidual(r);

		double beta = 0;
		Vector z = Precond.Solve(r);
		result.AddCost(Precond.SolvingComputationalWork());
		Vector d = z;
		this->IterationCount = 0;

		while (!StoppingCriteriaReached(result))
		{
			IterationResult result(result);

			if (this->IterationCount > 0)
				d = z + beta * d;

			double alpha = r.dot(z)/(d.dot(A*d));
			x = x + alpha * d;

			Vector old_r = r;
			Vector old_z = z;

			r = r - alpha * A * d;
			result.AddCost(2 * A.nonZeros());

			z = Precond.Solve(r);
			result.AddCost(Precond.SolvingComputationalWork());

			beta = (r.dot(z))/(old_r.dot(old_z));

			this->IterationCount++;

			result.SetX(x);
			result.SetResidual(r);

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		if (this->PrintIterationResults)
			cout << endl;

		return x;
	}
};