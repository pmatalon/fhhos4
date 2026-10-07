#pragma once
#include "IterativeSolver.h"
#include "../ProgramArguments.h"
using namespace std;

// U-AMG called through the library fhhos4_AMG (library/), as an external code calls it (fhhos4::Solver), to test the
// library against the U-AMG compiled in the program:
// - -s libuamg: U-AMG alone (as -s uamg);
// - -s fcglibuamg (or -s fcg -prec libuamg): the FCG of the library preconditioned by U-AMG (as -s fcguamg).
// The library does not give the cycle alone (U-AMG is a non-linear preconditioner, that requires FCG): libuamg cannot
// be the preconditioner of another solver of the program. The options of the program that fhhos4::Solver doesn't have
// keep their default values.
// Defined in LibraryUAMG.cpp, the only file of the program that includes the public header of the library: changing it
// doesn't recompile Program.cpp.
IterativeSolver* CreateLibraryUAMG(const ProgramArguments& args, int dim, int faceDegree, int cellDegree, bool useFCG);
