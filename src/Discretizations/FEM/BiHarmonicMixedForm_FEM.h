#pragma once
#include "../HHO/BiHarmonicMixedForm.h"
#include "Diffusion_FEM.h"
using namespace std;

template<int Dim>
class BiHarmonicMixedForm_FEM : public BiHarmonicMixedForm
{
public:
	BiHarmonicMixedForm_FEM() {}

	virtual Diffusion_FEM<Dim>& DiffPb() = 0;

	virtual Vector Solve1stDiffProblem(const Vector& bc) = 0;
	virtual Vector Solve2ndDiffProblem(const Vector& source, bool returnBoundaryOnly = false) = 0;

	virtual Vector Solve1stDiffProblem_Homogeneous(const Vector& bc) = 0;
	virtual Vector Solve2ndDiffProblem_Homogeneous(const Vector& source) = 0;

	Vector ProblemOperator(const Vector& x) override
	{
		Vector delta = Solve1stDiffProblem_Homogeneous(x);
		return -Solve2ndDiffProblem_Homogeneous(delta);
	}

	DenseMatrix Matrix()
	{
		Vector theta0 = FindCompatibleTheta();
		int n = theta0.rows();
		DenseMatrix A(n, n);

		int nThreads = 1;

		NumberParallelLoop<EmptyResultChunk> parallelLoop(n, nThreads);
		parallelLoop.Execute([this, n, &A](BigNumber i, ParallelChunk<EmptyResultChunk>* chunk)
			{
				Vector e_i = Vector::Zero(n);
				e_i[i] = 1;
				Vector lambda = Solve1stDiffProblem_Homogeneous(e_i);
				A.col(i) = -Solve2ndDiffProblem_Homogeneous(lambda);
			});
		return A;
	}
};
