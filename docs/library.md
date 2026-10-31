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
  (default) or Q, which gives the block sizes of the cells, faces and p-levels. `Result`, `Error`.
- `SetupSym()` from `A_TT`, `A_TF`, `A_FF` and the interpolations of 1, or `SetupNonSym()` with `A_FT` too,
  then `Solve(b_T, b_F, x_T, x_F)`: the library condenses the system by cells (`StaticCondensation.h`) and computes the
  cell unknowns. `CondensedMatrix()`, and `Solve(b, x)` on it (right-hand sides condensed by the caller).
- Matrices: `fhhos4::SparseMatrix` (Eigen row-major, of the library's index type `SparseIndex`, generated into
  `fhhos4_version.h` from `SMALL_INDEX`) or an Eigen map of one (`fhhos4::SparseMatrixMap`). The setups are templates on
  `Eigen::SparseCompressedBase<Derived>`: a wrong format is a `static_assert` naming the Eigen conversion
  `fhhos4::SparseMatrix B = A;`; a temporary `Eigen::SparseMatrix` is refused by a deleted overload where the solver
  keeps it (`A_FT`). Checked by one parallel pass, `CheckCSR()`.
- The solvers read the inputs as `BorrowedSparseMatrix` (`Utils/BorrowedSparseMatrix.h`: a `SparseMatrix` pointing
  Eigen's protected members to the caller's arrays). The library keeps `A_FT` (`SetupNonSym()`:
  `StaticCondensation` reads it at each solve); on several processes it computes its own from it (`docs/mpi.md`), but
  the rule for the caller stays the same.
- Helpers, documented as copies but one: `TransposedView()` (no copy: a symmetric column-major matrix as its transpose,
  an Eigen map); `ExtractSymmetricBlocks()`/`ExtractNonsymmetricBlocks()` (a global Eigen row-major matrix + `DoFLists`
  of the cells and faces, index types int32/int64, one parallel pass: `SparseMatrixOps::BuildByRows`; Eigen blocks,
  `HybridBlocks`, copied from the library's arrays); `Gather()`/`Scatter()`. `CondensedMatrix()`: an Eigen map. Other
  Eigen formats are converted by Eigen in the calling code.
- `fhhos4::detail` holds the implementation types, never in a signature the user calls nor in the docs:
  `BasicCSRView`/`CSRView` (the arrays given to the compiled functions), `CSRMatrix` (owned by the library), the
  compiled helpers (`detail::ExtractSymmetricBlocks()` on any index types, also mixed, `detail::ToCSR()`,
  `detail::Transpose()`), used by the Python module. The raw-view `Solver::Setup()` and `CondensedMatrixView()` are
  protected: the Python module's `PySolver` derives from `Solver`.
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
- Assemble into dedicated global matrices and vectors named as the blocks (`A_TT`, `A_TF`, `A_FT`, `A_FF`, `b_T`,
  `b_F`), never into the condensed system's variables: the assembly code must show that it assembles the blocks. Watch
  the memory: reserve, release, no restricted copies.
- The data prepared as if fhhos4 were the only solver; the solver part minimal: `SetupSym(A_TT, A_TF, A_FF,
  ...)` or `SetupNonSym(A_TT, A_TF, A_FT, A_FF, ...)`, then `Solve(b_T, b_F, x_T, x_F)`.
- No new parallel builders, adjacency lambdas or helper machinery in the client.
- The client's macros only around what does not compile without fhhos4: the include and the solver calls (with
  `#else`: an error, exit). The solver option and the block assembly stay unguarded.
- Check the client's bases: the interpolations of 1 (1/phi_0 for a hierarchical basis), and the order of the DoFs of
  each cell and face by degree (the p-levels extract the first ones).

## Threads, architecture, Eigen

- `Threads` default: at most one thread per physical core unless `OMP_NUM_THREADS` or OpenMP's binding is set (the
  hyperthreads slow the AMG down).
- Release builds use the conda compilers' `-march=nocona`.
- Eigen (ADR-0013): Eigen 3.4 only lacked `saveMarketDense` (`Utils/ExportModule.h` writes the format itself). Eigen 5's
  legacy macros say 3.5.0 (`EIGEN_WORLD_VERSION` stays 3; `EIGEN_PATCH_VERSION` exists from 5.0).
