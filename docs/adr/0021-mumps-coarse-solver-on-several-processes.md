# MUMPS for the coarsest level on several processes

On several processes the coarsest level of U-AMG is factorized and solved by MUMPS, distributed (2026-10-26, the
user: gathering the whole matrix on every process, each factorizing it, "is unacceptable if we want to achieve good
scaling on thousands of ranks"). The matrix, the right-hand side and the solution stay distributed: each process gives
its rows and gets back those of the solution.

- MUMPS is mandatory with `ENABLE_MPI` and optional without (`ENABLE_MUMPS`, off): ADR-0012 holds, since both options
  are off by default. MUMPS is CeCILL-C (LGPL-compatible) and the library links it privately (ADR-0001: still no GPL
  dependency). Its own dependencies depend on its build: conda-forge's `mumps-mpi` brings ParMETIS, whose license
  forbids commercial use; a MUMPS built without ParMETIS (with PT-Scotch) avoids it.
- `CoarseSolver::Cholesky` and `LU` keep their meaning: on several processes they are MUMPS's (`SYM=1`, `SYM=0`).
- On one process Eigen stays the default (the papers' results unchanged in every build); MUMPS is an explicit choice:
  `CoarseSolver::MUMPS`, `-coarse-solver mumpsch|mumpslu`.
