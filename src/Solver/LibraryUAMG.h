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
// The library is set up from the blocks of the uncondensed system (library/include/fhhos4.h): SetupFromBlocks(), then
// Solve(b_T, b_F, x_T, x_F): the library condenses the system and computes the cell unknowns. When the program needs the
// condensed system (see Program_Diffusion_HHO): Setup(A, blocks...), the same setup from the blocks (A ignored), then
// Solve(b) on the condensed matrix of the library.
// Defined in LibraryUAMG.cpp, the only file of the program that includes the public header of the library: changing it
// doesn't recompile Program.cpp.
class LibraryUAMG : public IterativeSolver
{
public:
	// Setup from the blocks of the uncondensed system (A_F_T empty if the matrix is symmetric: see
	// IterativeSolver::Setup()). The condensed matrix computed by the library becomes the matrix of the solver (Matrix).
	virtual void SetupFromBlocks(const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) = 0;

	// On several processes (the program's MPI_COMM_WORLD): from the rows of the cells and faces of the process (see
	// fhhos4::Solver::SetupSymmetric() and SetupNonsymmetric()): A_T_T, A_T_F (columns: the faces of the process, then its
	// ghosts, whose global numbers are faceGhosts), A_F_T (empty on all the processes if the matrix is symmetric; columns:
	// the cells of the process, then its ghost cells, cellGhosts), A_F_F (the rows of its faces, in full), the
	// interpolations of 1 on its cells and faces. Then Solve(b_T, b_F, x_T, x_F) with the vectors of the process.
	virtual void SetupFromDistributedBlocks(const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, const SparseMatrix& A_F_F, const vector<int64_t>& faceGhosts, const vector<int64_t>& cellGhosts, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) = 0;

	// Solves the uncondensed system, after SetupFromBlocks(). x_F: the initial guess on input.
	virtual void Solve(const Vector& b_T, const Vector& b_F, Vector& x_T, Vector& x_F) = 0;
	using IterativeSolver::Solve;
};

// krylov: the Krylov method of the library (fhhos4::Solver::Krylov); none: the cycles of U-AMG alone.
LibraryUAMG* CreateLibraryUAMG(const ProgramArguments& args, int dim, int faceDegree, int cellDegree, optional<KrylovMethod> krylov);
