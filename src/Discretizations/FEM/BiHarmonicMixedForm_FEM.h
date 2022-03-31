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
};
