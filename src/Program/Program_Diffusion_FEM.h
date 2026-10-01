#pragma once
#include "../ProgramArguments.h"

// Execute() is defined in Program_Diffusion_FEM_Impl.h, and compiled once in Program.cpp.
template <int Dim>
class Program_Diffusion_FEM
{
public:
	static void Execute(ProgramArguments& args);
};
