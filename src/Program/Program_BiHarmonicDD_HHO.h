#pragma once
#include "../ProgramArguments.h"

// Bi-harmonic equation in mixed form with Dirichlet BC enforced on both Laplacian problems (hence the suffix DD)
// Execute() is defined in Program_BiHarmonicDD_HHO_Impl.h, and compiled once in Program.cpp.
template <int Dim>
class Program_BiHarmonicDD_HHO
{
public:
	static void Execute(ProgramArguments& args);
};
