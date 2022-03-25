#pragma once
#include "BiHarmonicMixedForm.h"
#include "Diffusion_HHO.h"
#include "../TestCases/Diffusion/VirtualDiffusionTestCase.h"
#include "HigherOrderBoundary.h"
using namespace std;

template<int Dim>
class BiHarmonicMixedForm_HHO : public BiHarmonicMixedForm
{
public:
	BiHarmonicMixedForm_HHO() {}

	virtual Diffusion_HHO<Dim>& DiffPb() = 0;
};
