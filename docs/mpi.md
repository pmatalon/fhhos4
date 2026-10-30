# MPI: how it works in the code

Read before changing anything distributed (`ENABLE_MPI`, off by default). Decisions: `docs/adr/0011-mpi-design.md`.
Open work: `.scratch/mpi/`. Measurements: `PERFORMANCE.md` "Distributed memory (MPI)".

## Build and environment

MPI and MUMPS from the conda env `fhhos4` (`conda/environment.yml`: `mpich`, `mumps-mpi` and `_openmp_mutex=*=*_gnu`:
the pin keeps GNU's OpenMP, which MUMPS's build with OpenBLAS would replace by LLVM's; MUMPS then uses BLIS and
netlib's LAPACK: with several threads per process, set `BLIS_NUM_THREADS=1`, since BLIS takes `OMP_NUM_THREADS` and
creates its threads at each of MUMPS's small calls, see `PERFORMANCE.md`). Build directories: `build-mpi-lib`
(`-DENABLE_PROGRAM=OFF -DENABLE_MPI=ON`, ~1 min: the library and `fhhos4_mpi_tests`, the fast loop), `build-mpi` (the
program too: `-DENABLE_PROGRAM=ON`).

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
  Chebyshev, block Jacobi, (block) Gauss-Seidel = Gauss-Seidel in the process, Jacobi between them (`SweepRhs()`: not
  a contraction in the energy norm when the couplings between the processes are strong, e.g. 3D k=2, where FCG
  stagnates; the default smoother is Chebyshev, ADR-0017),
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
- Library: `Solver::Comm` (not set: one process, no MPI call; `MUMPS` then needs MPI initialized; forgotten on several
  processes: `ForgottenComm()` ends the errors of distributed inputs, ghosts or `A_FF` with more columns than rows), `SetupSymmetric(..., GhostNumbers faceGhosts)` (local columns: own faces then
  ghosts, one numbering for `A_TF` and `A_FF`, `A_FF` in full), `SetupNonsymmetric(..., faceGhosts, cellGhosts)`
  (`A_FT`: rows of the own faces, columns own cells then ghost cells; turned once by the setup into the rows of all
  the local faces, the columns of the own cells: `RowsWithGhosts`, kept by the condensation, which then runs as the
  symmetric one), `DistributedStaticCondensation.h`, `CollectiveCheck()` (errors of the inputs on all the processes). A
  process may own nothing (NGSolve 6.2.2608's rank 0). `SingularMatrix`: the kernel on the own faces, the dot products
  global, the pinned DoF of the coarsest level on one process (`docs/uamg.md`). Not yet: `CondensedMatrix()`. Python:
  `docs/python.md`.
- Program (`Program_Diffusion_HHO_Impl.h`, `Program/DistributedHybridSystem.h`): GMSH meshes only. The mesh of each
  process is its part (`GMSHMesh`, built when `Communicator::Current()` is distributed). The partition (`OpenPart()`:
  METIS, `Mesh.PartitionCreateGhostCells`, GMSH's split files `part_<rank + 1>.msh`) is kept in the cache directory
  (`-cache-dir`; shared by the nodes, unlike the temporary one), `<mesh cache name>_<P>parts`: read by the next runs on
  P processes (rank 0 then doesn't load the model); otherwise rank 0 loads the whole model (cache or generation),
  partitions it, writes it in `..._tmp<pid>` and renames it (atomic: no run reads an incomplete partition; if another
  run renamed first, its partition is kept), or deletes it after all have read their part (`-no-cache`). Collective
  checks: a partition with another number of parts, a part not found (the directory not shared). Each process opens
  its part, and keeps the ghost cells of GMSH (those sharing a node with a cell of the part) that share a face with one
  (`KeepGhostsOfFaces()`, `removeElements()` before `Build()`); in a part, the nodes are on entities of all dimensions,
  the entities of other parts empty. `ExchangeGhosts()` (after `Build()`, one request/reply to the processes of the
  ghosts, by GMSH tag): the global numbers of the cells (`ElementGlobalNumbers`, `ElementRanks`), the physical part of
  the ghosts (none in the files), and which of their faces with one cell in the part are on the boundary of the domain
  (with their boundary part): the others lie between the ghost and a cell of another part (interior faces with no
  `Element2`). Faces are identified between processes by the GMSH tags of their vertices (`FaceKey`), and oriented by
  them in a part (`OrientAlikeOnAllProcesses()`: the sign of the odd functions of the face bases, otherwise set by the
  first cell of the face, which differs between processes). `GlobalFaceNumbers()` (after the numbering of the faces
  by the problem, before the assembly): the faces of the system owned by the smaller rank of their cells, numbered
  by exclusive scan; the others asked to their process, first those of the cells of the process (both cells of such a
  face are in the part), then those of ghosts only, to the process of the ghost. `Diffusion_HHO::Restrict()`: the
  local matrices on all the cells of the part (the rows of its cells and faces complete, no exchange), the
  post-processing on its cells (`OwnElements()`), the sums reduced (L2 error, `HHOCellSpace`/`HHOReconstructSpace` and
  its measure, `HHOBoundarySpace::Integral()` on the boundary faces of its cells, `ZeroMeanEnforcer`,
  `NumericImageEnforcer`: the full Neumann right-hand side of the program's U-AMG projected on the rows of the
  process); `H()`, the regularity, the average h: of the whole mesh. `DistributedHybridSystem`: its rows in the global
  numbers; the solution of the faces of its rows (among which those of its cells) by a ghost update of `FaceMap`
  (`Local()`). `CheckDistributedArguments()` lists what is refused. `-uamg-nonsym`: `A_F_T` from the rows of
  the faces (`FromGlobalColumns()`: the ghost cells), local columns for the library, `RowsWithGhosts()` for the
  program's U-AMG.

## Tests

`fhhos4_mpi_tests` (`tests/MPI`, own `main()`, `ctest -R MPI_np`: 1-4 processes; `tests/support/PartitionedSystem.h`;
`MUMPSTest.cpp`: the coarse solver alone, with an empty first process, MUMPS's host);
`ctest -R ProgramMPI` (the program under `mpiexec`, L2 errors; 3D k=2 also the iterations of one process; `ProgramMPI_PartitionCache_*`: a partition written,
read by the next run, refused for another number of processes, in `build-mpi/test_partition_cache`); in an MPI build, `fhhos4_tests` uses
`tests/MPI/main.cpp` (MPI initialized), `ToySolver()` without `Comm`.
