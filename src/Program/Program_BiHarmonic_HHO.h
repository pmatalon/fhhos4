#pragma once
#include "../ProgramArguments.h"
#include "ProgramResults.h"

// Biharmonic equation in mixed form with mixed (homogeneous) Dirichlet-Neumann BC
// Execute() is defined in Program_BiHarmonic_HHO_Impl.h, and compiled once in Program.cpp.
template <int Dim>
class Program_BiHarmonic_HHO
{
public:
	static void Execute(ProgramArguments& args, ProgramResults* results = nullptr);
};
