# MPI: how it works in the code

Read before changing anything distributed (`ENABLE_MPI`, off by default). Decisions: `docs/adr/0011-mpi-design.md`.
Open work: `.scratch/mpi/`. Measurements: `PERFORMANCE.md` "Distributed memory (MPI)".

## Build and environment

MPI and MUMPS from the conda env (installed in `fhhos4`, not in `conda/environment.yml`: `mpich`, and
`conda install -c conda-forge "mumps-mpi=5.8.2" "_openmp_mutex=*=*_gnu"`: the pin keeps GNU's OpenMP, which MUMPS's
build with OpenBLAS would replace by LLVM's; MUMPS then uses BLIS and netlib's LAPACK: with several threads per
process, set `BLIS_NUM_THREADS=1`, since BLIS takes `OMP_NUM_THREADS` and creates its threads at each of MUMPS's small
calls, see `PERFORMANCE.md`). Build directories:
`build-mpi-lib` (`-DENABLE_PROGRAM=OFF -DENABLE_MPI=ON`, ~1 min: the library and `fhhos4_mpi_tests`, the fast loop),
`build-mpi` (the program too).

## Code

- `Utils/Communicator.h`: the MPI calls, the only file with `#ifdef MPI_ENABLED` but `Solver/LibraryUAMG.cpp` (the
  library's `Comm`; public macro: `FHHOS4_MPI` in `fhhos4_version.h`, `library/Solver.cpp`). `Communicator`
  (reductions, `AllToAll` for the setups, `Neighbours` for the halos, `FortranHandle()` for MUMPS),
  `Communicator::Current()` (the communicator of
  the reductions of the solvers: `MPI_COMM_WORLD` in the program, the solver's in a library call, `CallScope`),
  `Parallel::Dot/Norm/SquaredNorm` (on one process: Eigen's, bit for bit), `MPIEnvironment` (main(); `MPI_Finalize()`
  at exit only in the initializing process, not in GTest's death-test forks).
- `Parallelism::SharingProcesses()`: the processes of a node with the same affinity (mpiexec without binding) divide
  the default threads (OpenMP and solvers; the library: per `Comm`).
- On one process, nothing distributed is built: the sequential code runs (same results). On several:
  `Utils/IndexMap.h` (contiguous ownership in rank order, ghosts in any order, halo: `UpdateGhosts`,
  `AddGhostsToOwners`, `Forward/ReverseExchange`), `Utils/DistributedMatrix.h` (rows owned, `Diag` + `Offd` (ghost
  columns): the smoothers work on `Diag`; `MultiplyTransposed` for the restriction; `FetchGhostRows`, `RowsWithGhosts` (the rows of the ghosts fetched, the columns of another map),
  `AssembleFromLocalRows`: the rows of the ghosts added to their owners').
- Solvers: `Solver::Distributed`, `Setup(const DistributedMatrix&)` (FatalError unless overridden: FCG, GCR, Multigrid,
  Chebyshev, block Jacobi, (block) Gauss-Seidel = Gauss-Seidel in the process, Jacobi between them (`SweepRhs()`),
  Cholesky/LU by MUMPS (ADR-0021; `Solver/Direct/MUMPS.h`: each process gives its rows of the matrix and of the
  right-hand side; the solution, as MUMPS distributes it, goes to the owners of its rows by an `IndexMap` whose ghosts
  are the rows held for other processes, built at the first solve)); `IterativeSolver::SetDistributedMatrix()`; the blocks setup `Setup(A, A_T_T, A_T_F (local
  columns), A_F_T (empty if symmetric), faceMap, interps)`; `Level::DistributedOperator`, `DistributedP`. FCG applies the full matrix (not the lower
  triangle) on several processes.
- U-AMG: `UncondensedLevel::CoarsenInDistributedPasses()` with `DistributedPass.h`: each process aggregates its cells,
  remote cells/aggregates of the interface faces exchanged by the faces' halo, coarse faces between two processes owned
  by the smaller rank (numbered by request/reply), passes renumbered (owned, then ghosts), prolongations chained
  locally (interface faces are kept: injection rows known on both sides), operator assembled once per level. Only the
  default options (flat passes, chained prolongation, no `A_F_F`). Non-symmetric: `A_F_T` on the rows of all the local
  faces and the columns of the cells of the process (what its cells read: couplings, `P^T A_F_T Q_T`), the downwind
  order of the faces of each process (`Diag`; Jacobi between the processes). Ghost faces may belong to no local cell
  (couplings of `A_FF`): `CoarseningPass::NOwnedFaces`. All global decisions (coarsening ratio, errors) reduced.
- Library: `Solver::Comm` (no default), `SetupSymmetric(..., GhostNumbers faceGhosts)` (local columns: own faces then
  ghosts, one numbering for `A_TF` and `A_FF`, `A_FF` in full), `SetupNonsymmetric(..., faceGhosts, cellGhosts)`
  (`A_FT`: rows of the own faces, columns own cells then ghost cells; turned once by the setup into the rows of all
  the local faces, the columns of the own cells: `RowsWithGhosts`, kept by the condensation, which then runs as the
  symmetric one), `DistributedStaticCondensation.h`, `CollectiveCheck()` (errors of the inputs on all the processes). A
  process may own nothing (NGSolve 6.2.2608's rank 0). Not yet: `CondensedMatrix()`. Python: `docs/python.md`.
- Program (`Program_Diffusion_HHO_Impl.h`, `Program/DistributedHybridSystem.h`): GMSH meshes only, partitioned by GMSH
  (`GMSHMesh::PartitionElements()`, rank 0, broadcast); the mesh built by rank 0 first (cache); the system assembled by
  every process (replicated, not distributed yet), its rows extracted (faces owned by the smaller rank of their cells),
  the solution gathered; `CheckDistributedArguments()` lists what is refused. `-uamg-nonsym`: `A_F_T` from the rows of
  the faces (`FromGlobalColumns()`: the ghost cells), local columns for the library, `RowsWithGhosts()` for the
  program's U-AMG.

## Tests

`fhhos4_mpi_tests` (`tests/MPI`, own `main()`, `ctest -R MPI_np`: 1-4 processes; `tests/support/PartitionedSystem.h`;
`MUMPSTest.cpp`: the coarse solver alone, with an empty first process, MUMPS's host);
`ctest -R ProgramMPI` (the program under `mpiexec`, L2 errors); in an MPI build, `fhhos4_tests` uses
`tests/MPI/main.cpp` (MPI initialized), `ToySolver()` sets `Comm = MPI_COMM_SELF`.
