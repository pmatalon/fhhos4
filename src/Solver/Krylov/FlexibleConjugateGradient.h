#pragma once
#include "../IterativeSolver.h"
#include "../Preconditioner.h"
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
	SolverPreconditioner Precond;
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
	void Solve(const Vector& b, Vector& x, bool xEquals0, bool computeResidual, bool computeAx) override
	{
		const SparseMatrix& A = *this->Matrix;

		// Compared to CG, the work of 1 iteration of FCG is 1 more dot product by previous direction saved.
		// FCG(1) costs 1 more dot product by iteration than CG, i.e. 1 matrix-vector multiplication + 3 dot products.
		this->SolvingComputationalWork = 0; 

		if (this->ComputeExactSolution)
			this->ExactSolution = this->_directSolver.solve(b);

		IterationResult result = CreateFirstIterationResult(b, x);

		// Init restart parameters
		list<Direction> previousDirections;

		Vector& r = this->Residual;
		if (xEquals0)
		{
			r = b;
			result.SetResidualAsB();
		}
		else
		{
			r = b - A.selfadjointView<Eigen::Lower>() * x;            result.AddWorkInFlops(Cost::DAXPY(A));
			result.SetResidualNorm(r.norm());                         result.AddWorkInFlops(Cost::Norm(r));
		}

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			// Apply the preconditioner
			Vector z, Az;
			if (Precond.CanOptimizeResidualComputation())
			{
				std::tie(z, Az) = Precond.ApplyAndComputeAe(r);    result.AddWorkInMFlops(Precond.SolvingComputationalWork());
			}
			else
			{
				z = Precond.Apply(r);                              result.AddWorkInMFlops(Precond.SolvingComputationalWork());
			}

			// The final direction of research, d, is found by 
			// A-orthogonalizing z with the previous directions of research, 
			// according to the truncation / restart method chosen.

			Vector* d = new Vector(z);

			Vector* Ad = nullptr;
			if (Precond.CanOptimizeResidualComputation())
				Ad = new Vector(std::move(Az));

			for (Direction const& directionk : previousDirections)
			{
				Vector& dk  = *directionk.d;
				Vector& Adk = *directionk.Ad;
				double dk_dot_Adk = directionk.d_dot_Ad;

				// A-orthogonalization
				double z_dot_Adk = z.dot(Adk);                        result.AddWorkInFlops(Cost::Dot(z));
				*d -= (z_dot_Adk / dk_dot_Adk) * dk;                  result.AddWorkInFlops(Cost::VectorDAXPY(*d));

				if (Precond.CanOptimizeResidualComputation())
				{
					*Ad -= (z_dot_Adk / dk_dot_Adk) * Adk;            result.AddWorkInFlops(Cost::VectorDAXPY(*d));
				}
			}
			assert(d->norm() > 0);

			if (!Ad)
			{
				Ad = new Vector(A.selfadjointView<Eigen::Lower>() * (*d));  result.AddWorkInFlops(Cost::MatVec(A));
			}

			double d_dot_Ad = d->dot(*Ad);                            result.AddWorkInFlops(Cost::Dot(*d));
			assert(d_dot_Ad != 0);

			// Step length in the direction of research
			double alpha = r.dot(*d) / d_dot_Ad;                      result.AddWorkInFlops(Cost::Dot(r));

			// Moving from the current solution to the next
			// by taking the step in the direction of research
			x += alpha * (*d);                                        result.AddWorkInFlops(Cost::VectorDAXPY(x));

			// New residual
			r -= alpha * (*Ad);                                       result.AddWorkInFlops(Cost::VectorDAXPY(r));

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
			result.SetResidualNorm(r.norm());                         result.AddWorkInFlops(Cost::Norm(r));

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