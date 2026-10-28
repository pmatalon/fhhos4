# Distributed memory: the library first, GMSH partitions, aggregation per process

MPI (2026-10-20) targets U-AMG in the library, for DOLFINx (where MPI is mandatory); the program is its test bed.

- MPI is optional (`ENABLE_MPI`, off), and on one process the sequential code runs: nothing distributed is built.
- Only what U-AMG needs is distributed (Krylov methods, its smoothers, the coarse solver): not the geometric multigrid,
  the other AMGs, nor the biharmonic solvers.
- On several processes the program takes GMSH meshes only, partitioned by GMSH itself: no in-house partitioner, no
  METIS.
- Each process aggregates its own cells: the iteration counts may vary with the number of processes.
- The program assembles the matrix by recomputing rather than communicating: each process computes the local
  matrices of its own cells and of the other cells of its faces (a layer of ghost cells), so that the rows of its faces
  are complete with no exchange. The library, which receives rows and has no mesh, adds the rows of the ghosts to their
  owners.
- The library's layout is DOLFINx's, PETSc's and hypre's: each process gives the rows of its own cells and faces, the
  columns numbered locally (own, then ghosts), with the global numbers of its ghosts.
- `Solver::Comm` has no default (ADR-0005).
- Non-symmetric systems on several processes come after the symmetric case is complete.
