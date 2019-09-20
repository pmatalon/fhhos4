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
	}

private:
	Eigen::VectorXd Solve(const Eigen::VectorXd& b, Eigen::VectorXd& initialGuess) override
	{
		if (this->ComputeExactSolution)
			this->_exactSolution = this->_directSolver.solve(b);

		Eigen::VectorXd x = initialGuess;
		Eigen::VectorXd r = b - A * x;

		double beta = 0;
		Eigen::VectorXd z = Precond.Solve(r);
		Eigen::VectorXd d = z;
		this->IterationCount = 0;

		IterationResult result = SaveIterationResult(x, b, r);
		while (!StoppingCriteriaReached(result))
		{
			if (this->IterationCount > 0)
				d = z + beta * d;

			double alpha = r.dot(z)/(d.dot(A*d));
			x = x + alpha * d;

			Eigen::VectorXd old_r = r;
			Eigen::VectorXd old_z = z;
			r = r - alpha * A * d;
			z = Precond.Solve(r);

			beta = (r.dot(z))/(old_r.dot(old_z));

			this->IterationCount++;
			result = SaveIterationResult(x, b, r);
		}

		if (this->PrintIterationResults)
			cout << endl;

		return x;
	}
};