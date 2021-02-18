#pragma once
#include "IterativeSolver.h"
#include "Preconditioner.h"
#include <list>
using namespace std;

class FlexibleConjugateGradient : public IterativeSolver
{
private:
	struct Direction
	{
		Vector* d;
		Vector* Ad;
		double d_dot_Ad;
	};
public:
	Preconditioner Precond;
	int Truncation; // max number of previous search direction to use for A-orthogonalization
	bool Restart = false;

	FlexibleConjugateGradient(int truncation = 1, bool restart = false)
	{
		Truncation = truncation;
		Restart = restart;
	}

	virtual void Serialize(ostream& os) const override
	{
		os << "Flexible Conjugate Gradient, preconditioner: " << Precond;
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		this->Precond.Setup(A);
		this->SetupComputationalWork = this->Precond.SetupComputationalWork();
	}

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F) override
	{
		IterativeSolver::Setup(A);
		this->Precond.Setup(A, A_T_T, A_T_F, A_F_F);
		this->SetupComputationalWork = this->Precond.SetupComputationalWork();
	}

public:
	void Solve(const Vector& b, Vector& initialGuess, bool zeroInitialGuess) override
	{
		const SparseMatrix& A = *this->Matrix;

		// Compared to CG, the work of 1 iteration of FCG is 1 more dot product by previous direction saved.
		// FCG(1) costs 1 more dot product by iteration than CG, i.e. 1 matrix-vector multiplication + 3 dot products.
		this->SolvingComputationalWork = 0; 

		if (this->ComputeExactSolution)
			this->ExactSolution = this->_directSolver.solve(b);

		IterationResult result = CreateFirstIterationResult(b, initialGuess);

		// Init restart parameters
		list<Direction> previousDirections;

		Vector& x = initialGuess;
		Vector r;
		if (zeroInitialGuess)
			r = b;
		else
		{
			r = b - A * x;                                           result.AddCost(Cost::DAXPY(A));
		}
		result.SetResidual(r);

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			// Apply the preconditioner
			Vector z = Precond.Apply(r);                             result.AddCost(Precond.SolvingComputationalWork());

			// The final direction of research, d, is found by 
			// A-orthogonalizing z with the previous directions of research, 
			// according to the truncation / restart method chosen.

			Vector* d = new Vector(z);
			for (Direction const& directionk : previousDirections)
			{
				Vector& dk  = *directionk.d;
				Vector& Adk = *directionk.Ad;
				double dk_dot_Adk = directionk.d_dot_Ad;

				*d -= (z.dot(Adk) / dk_dot_Adk) * dk;                 result.AddCost(Cost::Dot(z));
			}
			assert(d->norm() > 0);

			Vector* Ad = new Vector(A * (*d));                        result.AddCost(Cost::MatVec(A));
			double d_dot_Ad = d->dot(*Ad);                            result.AddCost(Cost::Dot(*d));
			assert(d_dot_Ad != 0);

			// Step length in the direction of research
			double alpha = r.dot(*d) / d_dot_Ad;                      result.AddCost(Cost::Dot(r));

			// Moving from the current solution to the next
			// by taking the step in the direction of research
			x += alpha * (*d);

			// New residual
			r -= alpha * (*Ad);

			// Restart?
			if (previousDirections.size() == this->Truncation)
			{
				if (this->Restart)
					Clear(previousDirections);
				else
					DeleteOldest(previousDirections);
			}

			// The direction of research is saved 
			// (even after a restart, as advised by Notay)
			Direction direction = { d, Ad, d_dot_Ad };
			previousDirections.push_back(direction);

			//---------- End of iteration ----------//
			this->IterationCount++;

			result.SetX(x);
			result.SetResidual(r);

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		Clear(previousDirections);

		if (this->PrintIterationResults)
			cout << endl;

		this->SolvingComputationalWork = result.SolvingComputationalWork();
	}

private:
	void Clear(list<Direction>& previousDirections)
	{
		for (Direction& directionk : previousDirections)
		{
			delete directionk.d;
			delete directionk.Ad;
		}
		previousDirections.clear();
	}

	void DeleteOldest(list<Direction>& previousDirections)
	{
		Direction& oldestDirection = previousDirections.front();
		delete oldestDirection.d;
		delete oldestDirection.Ad;
		previousDirections.pop_front();
	}
};