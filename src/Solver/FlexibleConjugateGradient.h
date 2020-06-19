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
		Vector d;
		Vector Ad;
		double d_dot_Ad;
	};
	list<Direction> _previousDirections;
public:
	Preconditioner Precond;
	int Truncation; // max number of previous search direction to use for A-orthogonalization

	FlexibleConjugateGradient(int truncation = 1)
	{
		Truncation = truncation;
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

public:
	Vector Solve(const Vector& b, Vector& initialGuess) override
	{
		// Compared to CG, the work of 1 iteration of FCG is 1 more dot product by previous direction saved.
		// FCG(1) costs 1 more dot product by iteration than CG, i.e. 1 matrix-vector multiplication + 3 dot products.
		this->SolvingComputationalWork = 0; 

		if (this->ComputeExactSolution)
			this->_exactSolution = this->_directSolver.solve(b);

		IterationResult result = CreateFirstIterationResult(b, initialGuess);

		// Init restart parameters
		_previousDirections = list<Direction>();

		Vector x = initialGuess;
		Vector r = b - A * x;                                       result.AddCost(2 * A.nonZeros()); // Cost: 1 sparse MatVec
		result.SetResidual(r);

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			// Apply the preconditioner
			Vector z = Precond.Solve(r);                             result.AddCost(Precond.SolvingComputationalWork()); // Cost: preconditioner

			// The final direction of research, d, is found by 
			// A-orthogonalizing z with the previous directions of research, 
			// according to the truncation / restart method chosen.
			Vector d = z;
			for (Direction const& directionk : _previousDirections)
			{
				Vector dk         = directionk.d;
				Vector Adk        = directionk.Ad;
				double dk_dot_Adk = directionk.d_dot_Ad;

				d -= z.dot(Adk) / dk_dot_Adk * dk;                   result.AddCost(2 * z.rows());     // Cost: 1 Dot
			}

			Vector Ad = A * d;                                       result.AddCost(2 * A.nonZeros()); // Cost: 1 sparse MatVec
			double d_dot_Ad = d.dot(Ad);                             result.AddCost(2 * d.rows());     // Cost: 1 Dot

			// Step length in the direction of research
			double alpha = r.dot(d) / d_dot_Ad;                      result.AddCost(2 * r.rows());     // Cost: 1 Dot

			// Moving from the current solution to the next
			// by taking the step in the direction of research
			x = x + alpha * d;

			// New residual
			r = r - alpha * Ad;

			// Restart?
			if (_previousDirections.size() == this->Truncation)
				_previousDirections.clear();

			// The direction of research is saved 
			// (even after a restart, as advised by Notay)
			Direction direction = { d, Ad, d_dot_Ad };
			_previousDirections.push_back(direction);

			//---------- End of iteration ----------//
			this->IterationCount++;

			result.SetX(x);
			result.SetResidual(r);

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		_previousDirections.clear();

		if (this->PrintIterationResults)
			cout << endl;

		return x;
	}
};