# The library: how it is built

For the library's users: `README_library.md` and its pages, `library/doc/`. This file is for changing the library.
Decisions: `docs/adr/` (0002-0006, 0008, 0013, 0020, 0023). Its measurements: `PERFORMANCE.md`.

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

- `Solver`: parameters as public members; `Dimension`, `FaceDegree`, `CellDegree` (k - 1 to k + 1) and
  `PolynomialSpace` (P or Q), which give the block sizes of the cells, faces and p-levels: set, or inferred by the setup
  from the blocks (-1, `PolynomialSpace::Inferred`, the defaults: see below); `VectorComponents` multiplies them (a
  vector-valued unknown, interleaved: see `docs/uamg.md`; -1, the default: 1 when the three are set, refused
  otherwise); `CellBlockSize()`/`FaceBlockSize()` give these block sizes once set or inferred (`BlockSizes()`, the
  setup's computation). `Result`, `Error`.
- `SetupSym()` from `A_TT`, `A_TF`, `A_FF` and the near null space (below), or `SetupNonSym()` with `A_FT` too,
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

## The discretization inferred from the blocks

ADR-0023. `src/Solver/DiscretizationInference.h`, header-only, without MPI nor PETSc (also compiled into the PETSc
preconditioner); its messages start with `fhhos4: `, worded with the names of the members (`MemberNames`) or of
PETSc's options:
- `AnalyzeBlocks()`, on the blocks of a process: the DoFs per cell n_T, the first closed diagonal block of the stored
  entries of `A_TT` (full or lower-only), equal to the last one; the GCD of the boundaries of the runs of consecutive
  columns of the first row of each cell of `A_TF` (multiples of n_F; the cells in parallel with OpenMP, sequential in
  the PETSc preconditioner); whether a column is in the rows of two cells (a face between them, which makes that GCD
  n_F; usually among the first cells).
- `CombineAnalyses()`, the analyses of the processes: n_T the same on those with cells, the GCD of their GCDs, a face
  between two cells somewhere.
- `InferDiscretization()`: the readings (`Readings()`: the discretizations that the library takes with these DoFs and
  the parameters set), grouped by solver (`SolverKey()`: k >= 1, and the block sizes of the p-levels of degree >= 1
  below the finest); one solver resolves, several are an ambiguity, with the parameter to propose; `Message()`,
  `Description()`.

## The near null space

`cellNearNullSpace`, `faceNearNullSpace` (`detail::ConstMatrixView`: the columns of the caller's matrix, read in place),
after the block sizes (`Impl::ReduceNearNullSpace()`, collective: the same number of columns on all the processes,
so that they take the same path):
- As many columns as `VectorComponents`: they must be the constants of the components on hierarchical bases
  (`ConstantsOnly()`: the column c non-zero on the DoF c of each block only, to 1e-10 relative; proportional to the
  first column), reduced to the interpolation of 1 of each block, the first column's: what `UncondensedAMG` takes.
  Otherwise an error. The setups and the results are those of the former input, one coefficient per block.
- More columns (the rigid body modes): `Impl::NearNullSpaceBeyondConstants()`, `SetupSym()`, degrees >= 1.
  `LowDegreeFunctions()` keeps the rows of the DoFs of degree <= 1 of each cell and face (checked:
  the coefficients beyond the degree 1 zero to 1e-8 relative, the functions independent on each block; on several
  processes, those of the process, by `CollectiveCheck()`), given
  to the multigrid (`UncondensedAMG::CellNearNullSpace`, `FaceNearNullSpace`, `docs/uamg.md`) with ones as
  interpolations of 1. The blocks are read in place as in the other setups, and the whole solver works on the system
  as given: no change of basis outside the last p-level. `MixedPrecision` applies to the levels of degree >= 1 only
  (`docs/uamg.md`; `Verbosity` 2 prints it).
- `SingularMatrix`: the kernel is the functions given, those of the multigrid (`UncondensedAMG::Kernel()`, computed
  by its setup), which `Impl::CheckKernel()` checks after the setup of the multigrid and `SolveCondensed()` projects
  the right-hand side on.

## The inference in the setups

`library/Solver.cpp`, `InferDiscretization()`, in both setups after the CSR checks of the inputs and before
`CheckParameters()` (which needs the block sizes): the four parameters set, nothing; the three set, the space alone: P
if P and Q have the same DoFs, otherwise the analysis of each process, Q if it shows the DoFs of Q, P if not or if it
fails (the checks of the sizes and of the near null space then refuse a wrong one); otherwise the analysis of each
process (`CollectiveCheck()`), gathered (`AllGatherV()`), combined, and inferred. The values found are
written into the members (`DiscretizationParametersSet()` checks those set: `CellDegree` k - 1 to k + 1 on every path).

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
- Check the client's bases: the near null space is the client's own interpolation of 1 (hierarchical bases: one
  non-zero coefficient per block), and the order of the DoFs of each cell and face by degree (the p-levels extract the
  first ones).

## The examples

`examples/` (`library/doc/` refers to them; built with the tests, run by ctest): `fhhos4_solver.cpp`, the library
alone (`ExampleSolver`: on 3 processes with MPI; its exit code: converged), and those of PETSc (below). The C++ ones
share `HHO_diffusion.h`: an HHO discretization on n x n squares (degrees k and l, Legendre bases), the mesh distributed
by rows, but seen by the assemblies as unstructured (cells and faces by global numbers, the faces of each cell, the
cells of each face; HHO's meshes are general): no loop on rows. The Python one has its own copy. `fhhos4_solver.cpp`
assembles the distributed blocks itself (each process also computes the other cells of its faces, for the rows of its
faces): the library's own solve of `petsc_preconditioner.cpp` (`-ksp_type preonly -pc_fhhos4_krylov fcg`) gives its
iterations and residuals exactly.

## The preconditioner for PETSc

`library/petsc/fhhos4_petsc.cpp` (`ENABLE_PETSC_PC`, ADR-0012): `libfhhos4_petsc.so`, target `fhhos4_petsc` (exported
`fhhos4::petsc`), linked to `libfhhos4` (public) and to PETSc (private; pkg-config's `PETSc.pc`), installed with the
run path `$ORIGIN` (the programs that load it don't link `libfhhos4`: found next to it); its header
`library/include/fhhos4_petsc.h` (`PCFHHOS4`, `PetscDLLibraryRegister_fhhos4_petsc()`: the function PETSc calls when
`-dll_append` loads the library). One file in PETSc's style, on PETSc's interface for the implementations of
preconditioners (`petsc/private/pcimpl.h`), `fhhos4.h` and `src/Solver/DiscretizationInference.h` (header-only,
compiled into it: moving into PETSc, `src/ksp/pc/impls`, means installing this header too).
- Setup (`PCSetUp_FHHOS4()`): the operator a 2 x 2 MATNEST of AIJ blocks, the cells first, or an AIJ matrix with the
  index sets of its cells and faces composed on the PC (`"fhhos4_cells"`, `"fhhos4_faces"`: `PCFHHOS4SetIS()`), its
  blocks extracted by `MatCreateSubMatrix()` and destroyed at the end of the setup (the library keeps what it needs);
  from an index set that is not increasing, PETSc leaves the global numbers of a block's off-diagonal columns
  (`garray`) unsorted: `Block::SortOffd()` copies that part sorted, as `MergeOffDiag()` takes it. The rows in neither
  index set (`RowsLeftOut`): the identity, `VecCopy()` before the application. The options of the
  discretization not given inferred first (`InferDiscretization()`: `AnalyzeBlocks()` on the diagonal parts of `A_TT`
  and `A_TF`, its own `MPI_Allgather()`, the messages with the options' names; written into the solver's members; the
  dimension and the degrees given, the space alone, without analysis if P and Q have the same DoFs and solver, P if the
  stored entries don't give it, as `Solver::Setup()`); the
  block sizes from `Solver::CellBlockSize()`/`FaceBlockSize()`; the checks collective (`CheckOnAll()`: PETSc's errors
  are raised on the process that detects them, the others would wait in the next collective call): the DoFs of whole
  cells and faces on each process, the near null space, all its vectors on the DoFs of the process as the columns of the
  library's (`NearNullSpace()`; none: the constant one), which checks them. PETSc's arrays read in place (`MatGetRowIJ()`,
  `MatSeqAIJGetArrayRead()`; `MatMPIAIJGetSeqAIJ()` on several processes), released after the setup. One process: the
  blocks as Eigen maps of PETSc's arrays (copies by Eigen if `PetscInt` is not `SparseIndex`), `A_FT` copied (the solves
  read it). Several: `A_TT`'s diagonal part, the others by `MergeOffDiag()` (its errors local: reduced before the
  setup). The symmetric setup if `MatIsSymmetricKnown()` or `MatIsSPDKnown()`, or with `-pc_fhhos4_symmetric`.
- Application (`PCApply_FHHOS4()`): `Solve(r_T, r_F, z_T, z_F)` on the parts of the vectors (`VecGetSubVector()` with
  the index sets: plain vectors, or VECNEST ones, whose `VecGetSubVector()` takes the MATNEST's index sets; also for
  the near null space; with an AIJ operator's index sets, not contiguous, a copy through a scatter each time), `z_F`
  zero first; `Krylov::None`, `MaxIterations = 1` by default: one cycle.
- Options `-pc_fhhos4_*` (`PCSetFromOptions_FHHOS4()`), `PCView_FHHOS4()` (`-ksp_view`).
- Tests: `tests/PETSc/PCTest.cpp` (`fhhos4_petsc_tests`, its own `main()`: MPI, then PETSc; `ctest -R PETScPC`, 1-4
  processes, the library loaded by `-dll_append`, as by a user's program: not linked, it composes the index sets
  itself); the examples of `library/doc/petsc.md`,
  `examples/petsc_preconditioner.cpp` and `.py` (`ctest -R ExamplePETScPreconditioner`: `Cpp`, linked to the library,
  built with the tests; `Python`, if petsc4py is installed). Their HHO, the local matrix of the reference square, was
  checked by hand: symmetric, its kernel the interpolation of 1, exact on x^2; the C++ one (`examples/HHO_diffusion.h`,
  with the mesh and its distribution, shared with `examples/fhhos4_solver.cpp`) equals the Python one to rounding.
- Environment: conda env `fhhos4-petsc` (the stack of `fhhos4-python`, `docs/python.md`, with PETSc and petsc4py,
  pkg-config and GTest), build directory `build-python-petsc` (`-DENABLE_PROGRAM=OFF -DENABLE_PYTHON=ON -DENABLE_MPI=ON
  -DENABLE_PETSC_PC=ON -DENABLE_TESTS=ON`). On WSL, a Windows `pkg-config` on the PATH may be cached before the env's is
  installed: `cmake -UPKG_CONFIG_EXECUTABLE .`

## Threads, architecture, Eigen

- `Threads` default: at most one thread per physical core unless `OMP_NUM_THREADS` or OpenMP's binding is set (the
  hyperthreads slow the AMG down).
- Release builds use the conda compilers' `-march=nocona`.
- Eigen (ADR-0013): Eigen 3.4 only lacked `saveMarketDense` (`Utils/ExportModule.h` writes the format itself). Eigen 5's
  legacy macros say 3.5.0 (`EIGEN_WORLD_VERSION` stays 3; `EIGEN_PATCH_VERSION` exists from 5.0).
