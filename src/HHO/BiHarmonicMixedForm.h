#pragma once
#include "../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "../Solver/Solver.h"
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

	virtual Vector FindCompatibleTheta() = 0;

	virtual Vector Solve1stDiffProblemWithFSource(const Vector& bc) = 0;

	virtual Vector Solve1stDiffProblemWithZeroSource(const Vector& bc) = 0;

	virtual Vector Solve2ndDiffProblem(const Vector& source, bool returnBoundaryOnly = false) = 0;

	virtual double L2InnerProdOnBoundary(const Vector& v1, const Vector& v2) = 0;

	virtual Vector ComputeSolution(const Vector& theta) = 0;

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
	DenseMatrix Matrix()
	{
		Vector theta0 = FindCompatibleTheta();
		int n = theta0.rows();
		DenseMatrix A(n, n);
		DenseMatrix I = DenseMatrix::Identity(n, n);
		for (int i = 0; i < n; i++)
		{
			Vector lambda = Solve1stDiffProblemWithZeroSource(I.col(i));
			A.col(i) = -Solve2ndDiffProblem(lambda, true);
		}
		return A;
	}

	virtual ~BiHarmonicMixedForm() {}
};
