#pragma once
#include "ProgramArguments.h"
#include "Utils/Utils.h"
using namespace std;

// Prints the message as a fatal argument error and terminates the process.
// Shared by main()'s CLI parsing and by ApplyProgramArgumentDefaults() below.
inline void argument_error(string msg)
{
	cout << Utils::BeginRed << "Argument error: " << msg << Utils::EndColor << endl;
	cout << "------------------------- FAILURE -------------------------" << endl;
	exit(EXIT_FAILURE);
}

// Resolves every sentinel/"default" value left in a ProgramArguments (Dimension == -1,
// SolverCode/Mesher/MeshCode == "default", etc.) into a concrete, consistent configuration,
// exactly as main()'s CLI parsing does after getopt. Callers that build a ProgramArguments
// directly (e.g. tests) can call this instead of going through the CLI to get the same
// defaulting behavior.
//
// The boolean flags mirror whether the corresponding CLI option was explicitly passed by the
// user (main() tracks them while parsing getopt); left at their default of `true` (= "not
// explicitly set"), the same default-inference cascade as the CLI's own defaults applies.
void ApplyProgramArgumentDefaults(ProgramArguments& args,
	bool defaultRelativeCellPolyDegree = true,
	bool defaultTol2 = true,
	bool defaultCycle = true,
	bool defaultCoarseOperator = true,
	bool defaultCoarseSolver = true,
	bool defaultHPCoarseningStgy = true);
