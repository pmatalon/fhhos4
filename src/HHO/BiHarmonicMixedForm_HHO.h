#pragma once
#include "../TestCases/BiHarmonic/BiHarmonicTestCase.h"
#include "Diffusion_HHO.h"
#include "../TestCases/Diffusion/VirtualDiffusionTestCase.h"
#include "../Solver/Solver.h"
#include "HigherOrderBoundary.h"
using namespace std;

template<int Dim>
class BiHarmonicMixedForm_HHO
{
protected:
	Solver* _diffSolver = nullptr;

public:
	BiHarmonicMixedForm_HHO()
	{ }

	virtual Diffusion_HHO<Dim>& DiffPb() = 0;

	virtual void Setup() = 0;

	void SetDiffSolver(Solver* solver)
	{
		_diffSolver = solver;
		IterativeSolver* iter = dynamic_cast<IterativeSolver*>(_diffSolver);
		if (iter)
		{
			iter->PrintIterationResults = false;
			iter->MaxIterations = 50;
		}
	}

	virtual Vector FindCompatibleTheta() = 0;

	virtual Vector Solve1stDiffProblem(const Vector& neumann) = 0;

	virtual Vector Solve1stDiffProblemWithZeroSource(const Vector& neumann) = 0;

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
	void Matrix(ExportModule& out)
	{
		Vector theta0 = FindCompatibleTheta();
		int n = theta0.rows();
		DenseMatrix M(n, n);
		DenseMatrix I = DenseMatrix::Identity(n, n);
		for (int i = 0; i < n; i++)
		{
			Vector lambda = Solve1stDiffProblemWithZeroSource(I.col(i));
			M.col(i) = -Solve2ndDiffProblem(lambda, true);
		}
		//EigenSolver<MatrixXd> es;
		double det = M.determinant();
		auto v = M.eigenvalues();
		//cout << v << endl;
		out.ExportMatrix(M, "matrix");
	}

	virtual ~BiHarmonicMixedForm_HHO() {}
};
