# Distributed memory: the library first, GMSH partitions, aggregation per process

MPI (2026-10-20) targets U-AMG in the library, for DOLFINx (where MPI is mandatory); the program is its test bed.

- MPI is optional (`ENABLE_MPI`, off), and on one process the sequential code runs: nothing distributed is built.
- Only what U-AMG needs is distributed (Krylov methods, its smoothers, the coarse solver): not the geometric multigrid,
  the other AMGs, nor the biharmonic solvers.
- On several processes the program takes GMSH meshes only, partitioned by GMSH itself: no in-house partitioner, no
  METIS.
- Each process aggregates its own cells: the iteration counts may vary with the number of processes.
- The program assembles as DOLFINx (2026-11-07, the user): each process computes the local matrices of its own cells
  only; the rows of the faces of other processes, assembled locally, have their pattern sent once to the owners (whose
  rows get the couplings through the cells of other processes) and their values added in place to the owners' rows
  after each assembly. Its local numbering puts its own cells and faces first, so that the blocks it gives the library
  are the first rows of its matrices, not copies. Its mesh is its part, read from GMSH's partitioned files (one per
  process); GMSH's ghost cells only give the owners of the interface faces. "Recompute rather than communicate"
  (2026-10-28) was a rule of thumb, which standard practice overrides (the user, 2026-11-07): each process also
  computed the cells of other processes across its faces, so that the rows of its faces were complete with no
  exchange, at the price of cells computed twice and of what GMSH's files lack about these ghosts (their physical parts
  and boundary faces); no framework assembles that way. Until
  `.scratch/mpi/issues/12-petsc-binding.md` is done, the program still recomputes.
- The library's layout is the one of DOLFINx, PETSc and hypre after their assembly: each process gives the rows of its
  own cells and faces, complete, the columns numbered locally (own, then ghosts), with the global numbers of its
  ghosts; one CSR per block, where PETSc and hypre keep two (the owned and the ghost columns): their clients merge them,
  a copy. The library's condensation adds the contributions of the cells of each process to the rows of the faces of
  other processes to their owners.
- `Solver::Comm` not set: one process, with no MPI call (ADR-0005).
- Non-symmetric systems on several processes come after the symmetric case is complete.
