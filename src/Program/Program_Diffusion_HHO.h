#pragma once
#include "../ProgramArguments.h"
#include "ProgramResults.h"

// Execute() is defined in Program_Diffusion_HHO_Impl.h, and compiled once in Program.cpp.
template <int Dim>
class Program_Diffusion_HHO
{
public:
	static void Execute(ProgramArguments& args, ProgramResults* results = nullptr);
};
