#pragma once
#include <optional>
#include "IterativeSolver.h"
#include "../ProgramArguments.h"
using namespace std;

// U-AMG called through the library fhhos4 (library/), as an external code calls it (fhhos4::Solver), to test the
// library against the U-AMG compiled in the program:
// - -s libuamg: U-AMG alone (as -s uamg);
// - -s fcglibuamg (or -s fcg -prec libuamg): the FCG of the library preconditioned by U-AMG (as -s fcguamg);
// - -s gcrlibuamg (or -s gcr -prec libuamg): the GCR of the library preconditioned by U-AMG (as -s gcruamg).
// The library does not give the cycle alone (U-AMG is a non-linear preconditioner, that requires a flexible Krylov
// method): libuamg cannot be the preconditioner of another solver of the program. The options of the program that fhhos4::Solver doesn't have
// keep their default values.
// The two interfaces of the library (library/include/fhhos4.h):
// - the uncondensed system, by default: SetupFromBlocks() and Solve(b_T, b_F, x_T, x_F). The library condenses the
//   system and computes the cell unknowns;
// - the condensed system (-libuamg-condensed, or when the program needs it, see Program_Diffusion_HHO): the program
//   condenses the system, then Setup(A, blocks...) and Solve(b), as the codes that condense statically already.
// Defined in LibraryUAMG.cpp, the only file of the program that includes the public header of the library: changing it
// doesn't recompile Program.cpp.
class LibraryUAMG : public IterativeSolver
{
public:
	// Setup from the blocks of the uncondensed system (A_F_T empty if the matrix is symmetric: see
	// IterativeSolver::Setup()). The condensed matrix computed by the library becomes the matrix of the solver (Matrix).
	virtual void SetupFromBlocks(const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) = 0;

	// Solves the uncondensed system, after SetupFromBlocks(). x_F: the initial guess on input.
	virtual void Solve(const Vector& b_T, const Vector& b_F, Vector& x_T, Vector& x_F) = 0;
	using IterativeSolver::Solve;
};

// krylov: the Krylov method of the library (fhhos4::Solver::Krylov); none: the cycles of U-AMG alone.
LibraryUAMG* CreateLibraryUAMG(const ProgramArguments& args, int dim, int faceDegree, int cellDegree, optional<KrylovMethod> krylov);
