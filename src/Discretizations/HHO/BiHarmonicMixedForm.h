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

	virtual Vector Solve1stDiffProblem(const Vector& bc) = 0;
	virtual Vector Solve2ndDiffProblem(const Vector& source, bool returnBoundaryOnly = false) = 0;

	virtual Vector Solve1stDiffProblem_Homogeneous(const Vector& bc) = 0;
	virtual Vector Solve2ndDiffProblem_Homogeneous(const Vector& source) = 0;

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
	virtual DenseMatrix Matrix()
	{
		Vector theta0 = FindCompatibleTheta();
		int n = theta0.rows();
		DenseMatrix A(n, n);
		Vector e_i = Vector::Zero(n);

		NumberParallelLoop<EmptyResultChunk> parallelLoop(n);
		parallelLoop.Execute([this, n, &A](BigNumber i, ParallelChunk<EmptyResultChunk>* chunk)
			{
				Vector e_i = Vector::Zero(n);
				e_i[i] = 1;
				Vector lambda = Solve1stDiffProblem_Homogeneous(e_i);
				A.col(i) = -Solve2ndDiffProblem_Homogeneous(lambda);
			});
		return A;
	}

	virtual ~BiHarmonicMixedForm() {}
};
