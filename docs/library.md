# The library: how it is built

For the library's users: `README_library.md`. This file is for changing the library. Decisions: `docs/adr/` (0002-0006,
0008, 0013, 0020). Its measurements: `PERFORMANCE.md`.

## Target and files

- Shared library `libfhhos4.so`, CMake target `fhhos4_lib` (`fhhos4` is the program), exported as `fhhos4::fhhos4`.
  Compiled with `-fvisibility=hidden`: only the public API is exported.
- Public interface `library/include/fhhos4.h` (`#include <fhhos4>`, a file that includes `fhhos4.h`, or `<fhhos4.h>`;
  namespace `fhhos4`, pimpl); `Error` in `fhhos4_error.h`, which `Utils.h` includes alone. The public headers start
  with an `SPDX-License-Identifier` line and their copyright line.
- Implementation `library/Solver.cpp`. It builds the solver as the program does (`AlgebraicSolverFactory`,
  `ApplyUncondensedAMGDefaults`): its path must not include mesh or discretization headers, nor read the global
  `Utils::ProgramArgs`.

## Interface

- `Solver`: parameters as public members; `Dimension`, `FaceDegree` and `CellDegree` must be set; `PolynomialSpace`: P
  (default) or Q, which gives the block sizes of the cells, faces and p-levels; `VectorComponents` (default 1)
  multiplies them (a vector-valued unknown, interleaved: see `docs/uamg.md`); `CellBlockSize()`/`FaceBlockSize()`
  give these block sizes (`BlockSizes()`, the setup's computation). `Result`, `Error`.
- `SetupSym()` from `A_TT`, `A_TF`, `A_FF` and the interpolations of 1, or `SetupNonSym()` with `A_FT` too,
  then `Solve(b_T, b_F, x_T, x_F)`: the library condenses the system by cells (`StaticCondensation.h`) and computes the
  cell unknowns. `CondensedMatrix()` (one process), and `Solve(b, x)` on it (right-hand sides condensed by the caller).
  `CondenseRHS()`, `ApplyCondensedMatrix()` (`y = A x`: `Impl::A`, or `DistributedA` on several processes) and
  `SolveCellUnknowns()` (the steps of `Solve(b_T, b_F, x_T, x_F)`, `StaticCondensation`): for a Krylov method of the
  caller on the condensed system, preconditioned by one cycle (`Krylov::None`, `MaxIterations = 1`), on any number of
  processes.
- Matrices: `fhhos4::SparseMatrix` (Eigen row-major, of the library's index type `SparseIndex`, generated into
  `fhhos4_version.h` from `ENABLE_64BIT_INDEX`) or an Eigen map of one (`fhhos4::SparseMatrixMap`). The setups are
  templates on `Eigen::SparseCompressedBase<Derived>`: a wrong format is a `static_assert` naming the Eigen conversion
  `fhhos4::SparseMatrix B = A;`; a temporary `Eigen::SparseMatrix` is refused by a deleted overload where the solver
  keeps it (`A_FT`). Checked by one parallel pass, `CheckCSR()`; `A_TT` block diagonal (`CheckBlockDiagonal()`), and
  the triangles the setup reads present (`CheckTriangles()`: a lower-only `A_FF` where it is read in full, on several
  processes or by `SetupNonSym()`, would silently give another matrix), one parallel pass each.
- The solvers read the inputs as `BorrowedSparseMatrix` (`Utils/BorrowedSparseMatrix.h`: a `SparseMatrix` pointing
  Eigen's protected members to the caller's arrays). The library keeps `A_FT` (`SetupNonSym()`:
  `StaticCondensation` reads it at each solve); on several processes it computes its own from it (`docs/mpi.md`), but
  the rule for the caller stays the same.
- Helpers, documented as copies but one: `TransposedView()` (no copy: a symmetric column-major matrix as its transpose,
  an Eigen map); `ExtractBlocksSym()`/`ExtractBlocksNonSym()` (a global Eigen row-major matrix + `DoFLists`
  of the cells and faces, index types int32/int64, one parallel pass: `SparseMatrixOps::BuildByRows`; Eigen blocks,
  `HybridBlocks`, copied from the library's arrays); `SplitVector()`/`JoinVector()`; `MergeOffDiag()` (on
  several processes, the blocks in two parts, diagonal and off-diagonal, as PETSc's MPIAIJ and hypre's ParCSR: one
  parallel pass, no communication, no sort; the ghosts, the union of those of the blocks' off-diagonal columns;
  `DistributedBlocks`, ghosts as `std::vector`; `A_TT` not merged, its diagonal part read in place). `CondensedMatrix()`:
  an Eigen map. Other Eigen formats are converted by Eigen in the calling code.
- `fhhos4::detail` holds the implementation types, never in a signature the user calls nor in the docs:
  `BasicCSRView`/`CSRView` (the arrays given to the compiled functions; `LibraryView()`, an Eigen matrix in the
  format of the setups), `CSRMatrix` and `IndexArray` (owned by the library), `TwoPartBlock` (the parts of a block of
  `MergeOffDiag()`, one index type for both, int32/int64, and the global numbers, int32/int64), the compiled helpers
  (`detail::ExtractBlocksSym()` on any index types, also mixed, `detail::MergeOffDiag()`, `detail::ToCSR()`,
  `detail::Transpose()`), used by the Python module. The raw-view `Solver::Setup()` and
  `CondensedMatrixView()` are protected: the Python module's `PySolver` derives from `Solver`.
- The compiled functions take plain types only (views, `detail::` vector views); the Eigen functions are inline in the
  header, compiled with the client's Eigen.
- Parameters with a choice are enums (`fhhos4::Krylov`, `Cycle`, `PCoarsening`, `Smoother`, `Prolongation`,
  `CoarseSolver`; no strings), as the program's multigrid parameters (`SmootherType`, `CycleType`, `P_CoarsStgy`,
  `UAMGProlongation` in `Utils/Types.h`; the CLI codes in `SmootherFactory`, the coarse solver stays a solver code):
  `Solver/LibraryParameters.h` maps them both ways (one table per parameter; `CoarseSolver::MUMPS`: `mumpsch` or
  `mumpslu` for the symmetry of the matrix), for `library/Solver.cpp` and `LibraryUAMG.cpp`.

## In the program

`-s libuamg|fcglibuamg|gcrlibuamg` (`Solver/LibraryUAMG.cpp`, the only file of the program that includes the public
header, through `LibraryParameters.h`) runs the library on the uncondensed system (`Diffusion_HHO::CondensationBySolver`:
no static condensation in the program). For the exports of the system or of the algebraic error, the program
condenses the system, sets the library up from the blocks and solves the condensed system by `Solve(b, x)`
(`LibraryUAMG::Setup()`; tested through `-export lsys`). The full Neumann problem: `SingularMatrix`
(`SolverFactory::IsFullNeumann()`, also for `-s uamg`), the library projects b (the program does when it condenses);
on several processes too (see `docs/uamg.md`). `tests/UAMGLibraryTest.cpp` compares it with
`uamg|fcguamg|gcruamg`. `-uamg-nonsym` runs the non-symmetric path on the symmetric matrices; `-uamg-mixed` and
`-uamg-bcsr` the mixed precision and the block CSR storage (off in the program: the papers).

## Integrating the library into a client code

- Reuse what the client already computes and where it assembles: its local matrices, its assembly loops (the blocks
  assembled there, beside its own matrices), its Dirichlet treatment (where it moves the boundary values to the
  right-hand side).
- Assemble into dedicated matrices and vectors named as the blocks (`A_TT`, `A_TF`, `A_FT`, `A_FF`, `b_T`,
  `b_F`), never into the condensed system's variables: the assembly code must show that it assembles the blocks. Watch
  the memory: reserve, release, no restricted copies.
- The data prepared as if fhhos4 were the only solver; the solver part minimal: `SetupSym(A_TT, A_TF, A_FF,
  ...)` or `SetupNonSym(A_TT, A_TF, A_FT, A_FF, ...)`, then `Solve(b_T, b_F, x_T, x_F)`.
- No new parallel builders, adjacency lambdas or helper machinery in the client.
- The client's macros only around what does not compile without fhhos4: the include and the solver calls (with
  `#else`: an error, exit). The solver option and the block assembly stay unguarded.
- Check the client's bases: the interpolations of 1 (1/phi_0 for a hierarchical basis), and the order of the DoFs of
  each cell and face by degree (the p-levels extract the first ones).

## The examples

`examples/` (`README_library.md` refers to them; built with the tests, run by ctest): `fhhos4_solver.cpp`, the library
alone (`ExampleSolver`: on 3 processes with MPI; its exit code: converged), and those of PETSc (below). The C++ ones
share `HHO.h`: an HHO discretization on n x n squares (degrees k and l, Legendre bases), the mesh distributed by rows.
`fhhos4_solver.cpp` assembles the distributed blocks itself (each process also computes the cells of the row above its
own, for the rows of its faces): the library's own solve of `petsc_preconditioner.cpp` (`-ksp_type preonly
-pc_fhhos4_krylov fcg`) gives its iterations and residuals exactly.

## The preconditioner for PETSc

`library/petsc/fhhos4_petsc.cpp` (`ENABLE_PETSC_PC`, ADR-0012): `libfhhos4_petsc.so`, target `fhhos4_petsc` (exported
`fhhos4::petsc`), linked to `libfhhos4` (public) and to PETSc (private; pkg-config's `PETSc.pc`), installed with the
run path `$ORIGIN` (the programs that load it don't link `libfhhos4`: found next to it); its header
`library/include/fhhos4_petsc.h` (`PCFHHOS4`, `PetscDLLibraryRegister_fhhos4_petsc()`: the function PETSc calls when
`-dll_append` loads the library). One file in PETSc's style, on PETSc's interface for the implementations of
preconditioners (`petsc/private/pcimpl.h`) and `fhhos4.h` only, to move into PETSc (`src/ksp/pc/impls`).
- Setup (`PCSetUp_FHHOS4()`): the operator a 2 x 2 MATNEST of AIJ blocks, the cells first; the block sizes from
  `Solver::CellBlockSize()`/`FaceBlockSize()`; the checks collective (`CheckOnAll()`: PETSc's errors are raised on the
  process that detects them, the others would wait in the next collective call): the DoFs of whole cells and faces on
  each process, the interpolations of 1 from the near null space (its first vector, or the constant one: the first DoF
  of each block; the others 0, to 1e-10 relative). PETSc's arrays read in place (`MatGetRowIJ()`,
  `MatSeqAIJGetArrayRead()`; `MatMPIAIJGetSeqAIJ()` on several processes), released after the setup. One process: the
  blocks as Eigen maps of PETSc's arrays (copies by Eigen if `PetscInt` is not `SparseIndex`), `A_FT` copied (the solves
  read it). Several: `A_TT`'s diagonal part, the others by `MergeOffDiag()` (its errors local: reduced before the
  setup). The symmetric setup if `MatIsSymmetricKnown()` or `MatIsSPDKnown()`, or with `-pc_fhhos4_symmetric`.
- Application (`PCApply_FHHOS4()`): `Solve(r_T, r_F, z_T, z_F)` on the parts of the vectors (`VecGetSubVector()` with
  the nest's index sets: plain vectors, or VECNEST ones, whose `VecGetSubVector()` takes the MATNEST's index sets; also
  for the near null space), `z_F` zero first; `Krylov::None`, `MaxIterations = 1` by default: one cycle.
- Options `-pc_fhhos4_*` (`PCSetFromOptions_FHHOS4()`), `PCView_FHHOS4()` (`-ksp_view`).
- Tests: `tests/PETSc/PCTest.cpp` (`fhhos4_petsc_tests`, its own `main()`: MPI, then PETSc; `ctest -R PETScPC`, 1-4
  processes, the library loaded by `-dll_append`, as by a user's program); the examples of `README_library.md`,
  `examples/petsc_preconditioner.cpp` and `.py` (`ctest -R ExamplePETScPreconditioner`: `Cpp`, linked to the library,
  built with the tests; `Python`, if petsc4py is installed). Their HHO, the local matrix of the reference square, was
  checked by hand: symmetric, its kernel the interpolation of 1, exact on x^2; the C++ one (`examples/HHO.h`, with the
  mesh and its distribution, shared with `examples/fhhos4_solver.cpp`) equals the Python one to rounding.
- Environment: conda env `fhhos4-petsc` (the stack of `fhhos4-python`, `docs/python.md`, with PETSc and petsc4py,
  pkg-config and GTest), build directory `build-python-petsc` (`-DENABLE_PROGRAM=OFF -DENABLE_PYTHON=ON -DENABLE_MPI=ON
  -DENABLE_PETSC_PC=ON`). On WSL, a Windows `pkg-config` on the PATH may be cached before the env's is installed:
  `cmake -UPKG_CONFIG_EXECUTABLE .`

## Threads, architecture, Eigen

- `Threads` default: at most one thread per physical core unless `OMP_NUM_THREADS` or OpenMP's binding is set (the
  hyperthreads slow the AMG down).
- Release builds use the conda compilers' `-march=nocona`.
- Eigen (ADR-0013): Eigen 3.4 only lacked `saveMarketDense` (`Utils/ExportModule.h` writes the format itself). Eigen 5's
  legacy macros say 3.5.0 (`EIGEN_WORLD_VERSION` stays 3; `EIGEN_PATCH_VERSION` exists from 5.0).
