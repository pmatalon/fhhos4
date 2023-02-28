#pragma once
#include "BiHarmonicMixedForm.h"
#include "Diffusion_HHO.h"
using namespace std;

template<int Dim>
class BiHarmonicMixedForm_HHO : public BiHarmonicMixedForm
{
public:
	BiHarmonicMixedForm_HHO() {}

	virtual Diffusion_HHO<Dim>& DiffPb() = 0;

	virtual Vector Solve1stDiffProblem(const Vector& bc) = 0;
	virtual Vector Solve2ndDiffProblem(const Vector& source, bool returnBoundaryOnly = false) = 0;

	virtual DenseMatrix Matrix() = 0;
};
