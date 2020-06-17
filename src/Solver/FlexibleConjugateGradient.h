#pragma once
#include "IterativeSolver.h"
#include "Preconditioner.h"
#include <list>
using namespace std;

class FlexibleConjugateGradient : public IterativeSolver
{
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

private:
	//int _remainingIterBeforeRestart;
	list<Vector> _directions;
	//int _nStoredDirections; // truncation number

	Vector Solve(const Vector& b, Vector& initialGuess) override
	{
		this->SolvingComputationalWork = 0;

		if (this->ComputeExactSolution)
			this->_exactSolution = this->_directSolver.solve(b);

		IterationResult result = CreateFirstIterationResult(b, initialGuess);

		// Init restart parameters
		_directions = list<Vector>();

		Vector x = initialGuess;
		Vector r = b - A * x;
		result.AddCost(2 * A.nonZeros());
		result.SetResidual(r);

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			cout << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			// Apply the preconditioner
			Vector z = Precond.Solve(r);
			result.AddCost(Precond.SolvingComputationalWork());

			// The final direction of research is found by A-orthogonalizing z with the previous directions of research, 
			// according to the truncation / restart method chosen.
			Vector d = AOrthogonalizeToPreviousDirections(z);
			result.AddCost(_directions.size()*(2 * A.nonZeros()));

			Vector Ad = A * d;
			result.AddCost(2 * A.nonZeros());

			// Step length in the direction of research
			double alpha = r.dot(d)/(d.dot(Ad));

			// Moving from the current solution to the next by taking the step in the direction of research
			x = x + alpha * d;

			// New residual
			r = r - alpha * Ad;

			// Restart?
			if (_directions.size() == this->Truncation)
				_directions.clear();

			// The direction of research is saved (even after a restart, as advised by Notay)
			_directions.push_back(d);

			//---------- End of iteration ----------//
			this->IterationCount++;

			result.SetX(x);
			result.SetResidual(r);

			if (this->PrintIterationResults)
				cout << result << endl;
		}

		_directions.clear();

		if (this->PrintIterationResults)
			cout << endl;

		return x;
	}

	inline Vector AOrthogonalizeToPreviousDirections(const Vector& d)
	{
		return d - ProjectOnPreviousDirections(d);
	}

	Vector ProjectOnPreviousDirections(const Vector& z)
	{
		Vector result = Vector::Zero(z.rows());
		for (Vector const& dk : _directions)
		{
			Vector Adk = this->A * dk;
			result = result + (z.dot(Adk))/(dk.dot(Adk)) * dk;
		}
		return result;
	}
};