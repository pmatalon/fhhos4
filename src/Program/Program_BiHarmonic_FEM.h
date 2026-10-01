#pragma once
#include "../ProgramArguments.h"

// Biharmonic equation in mixed form with mixed (homogeneous) Dirichlet-Neumann BC
// Execute() is defined in Program_BiHarmonic_FEM_Impl.h, and compiled once in Program.cpp.
template <int Dim>
class Program_BiHarmonic_FEM
{
public:
	static void Execute(ProgramArguments& args);
};
