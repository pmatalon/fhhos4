#pragma once
#include "ProgramArguments.h"
using namespace std;

class Program
{
public:
	Program() {}
	virtual void Start(ProgramArguments& args) = 0;
	virtual ~Program() {}
};

// The members are defined in Program.cpp and explicitly instantiated there for the enabled
// dimensions, so that the heavy headers they need are compiled only once.
template <int Dim>
class ProgramDim : public Program
{
public:
	ProgramDim() : Program() {}

	void Start(ProgramArguments& args) override;

	// Sets the global state read by the program (Utils::ProgramArgs, mesh directories, GMSH
	// settings) from the arguments. Called by Start().
	static void InitGlobalState(const ProgramArguments& args);
};
