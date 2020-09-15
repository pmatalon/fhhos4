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
		const SparseMatrix& A = *this->Matrix;

		this->SolvingComputationalWork = 0; // The work of one iteration of CG is 1 matrix-vector product + 2 dot products

		if (this->ComputeExactSolution)
			this->_exactSolution = this->_directSolver.solve(b);

		IterationResult result = CreateFirstIterationResult(b, initialGuess);

		Vector x = initialGuess;
		Vector r = b - A * x;                               result.AddCost(2 * A.nonZeros()); // Cost: 1 sparse MatVec
		result.SetResidual(r);

		double beta = 0;
		Vector z = Precond.Solve(r);                        result.AddCost(Precond.SolvingComputationalWork()); // Cost: preconditioner
		Vector d = z;
		double r_dot_z = r.dot(z);                          result.AddCost(2 * r.rows());    // Cost: 1 Dot
		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			if (this->IterationCount > 0)
				d = z + beta * d;

			Vector Ad = A * d;                                  result.AddCost(2 * A.nonZeros()); // Cost: 1 sparse MatVec
			double alpha = r_dot_z / (d.dot(Ad));               result.AddCost(2 * r.rows());     // Cost: 1 Dot
			x = x + alpha * d;

			double old_r_dot_old_z = r_dot_z; // save the dot product before overwriting r and z

			r = r - alpha * Ad;
			z = Precond.Solve(r);                               result.AddCost(Precond.SolvingComputationalWork()); // Cost: preconditioner

			r_dot_z = r.dot(z);                                 result.AddCost(2 * r.rows());     // Cost: 1 Dot
			beta = r_dot_z / old_r_dot_old_z;

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