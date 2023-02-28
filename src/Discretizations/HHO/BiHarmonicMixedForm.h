#pragma once
#include "../../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "../../Solver/IterativeSolver.h"
using namespace std;

class BiHarmonicMixedForm
{
protected:
	Solver* _diffSolver = nullptr;

public:
	BiHarmonicMixedForm() {}

	virtual void Setup() = 0;

	void SetDiffSolver(Solver* solver)
	{
		_diffSolver = solver;
		IterativeSolver* iter = dynamic_cast<IterativeSolver*>(_diffSolver);
		if (iter)
		{
			iter->PrintIterationResults = false;
			iter->StoppingCrit = StoppingCriteria::NormalizedResidual;
			iter->MaxIterations = 50;
		}
	}

	IterativeSolver* IterativeDiffSolver()
	{
		return dynamic_cast<IterativeSolver*>(_diffSolver);
	}

	void SetDiffSolverTolerance(double tol)
	{
		IterativeSolver* iter = dynamic_cast<IterativeSolver*>(_diffSolver);
		if (iter)
			iter->Tolerance = tol;
	}

	double DiffSolverTolerance()
	{
		IterativeSolver* iter = dynamic_cast<IterativeSolver*>(_diffSolver);
		return iter ? iter->Tolerance : 0;
	}

	virtual Vector FindCompatibleTheta() = 0;

	virtual Vector ProblemOperator(const Vector& bc) = 0;

	virtual double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2) = 0;

	virtual pair<Vector, Vector> ComputeSolution(const Vector& theta) = 0;

protected:
	void CheckDiffSolverConvergence()
	{
		IterativeSolver* iterSolver = dynamic_cast<IterativeSolver*>(_diffSolver);
		if (iterSolver)
		{
			if (iterSolver->IterationCount == iterSolver->MaxIterations)
				Utils::Warning("The diffusion solver has reached the max number of iterations (" + to_string(iterSolver->MaxIterations) + ")");
		}
	}

public:

	virtual ~BiHarmonicMixedForm() {}
};
