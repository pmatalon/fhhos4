# Reproducibility modes: one option for the trade-off between speed and reproducible solutions

One option (2026-11-26, the user), `fhhos4::Solver::Reproducibility` (program `-reproducibility`, PETSc
`-pc_fhhos4_reproducibility`), chooses the reproducibility that the solver guarantees; the parallel algorithms that
would break it are replaced, the fastest kept within it. The guarantee: the same solution, bit for bit, from run to run,
given the same binary, machine, matrix and parameters, and the same number of processes (each process aggregates its
own cells, ADR-0011: the solution depends on the partition).

- `None` (the default, in the library and in the program): no guarantee. On several processes MUMPS gives solutions
  that differ at the rounding level from run to run (its scheduling between the processes; `ICNTL(28)` changes
  nothing, measured 2026-11-25); today's behaviour, unchanged.
- `OnFixedNThreads`: for the same number of threads per process. On several processes, the coarse matrix is
  gathered on the process 0 and factorized there by Eigen (`CoarseSolver` `Default`, `Cholesky`, `LU`), its
  right-hand sides gathered and its solutions scattered; `CoarseSolver::MUMPS` is refused there (on one process it is
  reproducible, and allowed).
- `OnAnyNThreads`: whatever the number of threads. Also the default smoother for a non-symmetric matrix becomes the
  l1 block Gauss-Seidel (`L1BlockGaussSeidel`, `l1bgs`: sequential in each process, l1 on the rows coupled to other
  processes), and the hybrid smoothers set explicitly are refused (their results depend on the number of threads,
  ADR-0018, and so does their 20% switch to sequential sweeps).

Chebyshev (ADR-0017), the sequential Gauss-Seidel sweeps, the setup and the Krylov methods (dot products by blocks of
rows in a fixed order) already give the same results on any number of threads: the strict modes cost little today
(the non-symmetric sweeps on one thread per process, the serial coarse solve on several processes, `PERFORMANCE.md`).
The option frames the future parallel algorithms: a split greedy aggregation among the threads would be allowed in
`None` only, or also in `OnFixedNThreads`.

The gathered coarse solve against ADR-0021: the process 0 holds and factorizes the whole coarse matrix, a serial step
that doesn't scale to thousands of processes, the reason ADR-0021 chose MUMPS. Acceptable because it is opt-in: ADR-0021
stays true for `None`. The matrix is gathered in rank order, its rows in global numbering, so that it doesn't depend
on the order in which the messages arrive.

Not covered, documented: MPI's reductions (`MPI_Allreduce`), deterministic in practice but not by the MPI standard;
the program's OpenMP `reduction(+)` on doubles outside the solve (the L2 errors, the full Neumann compatibility check).

## Considered Options

- "Up to rounding" modes (the user's first list, five modes): a rounding difference can change an iteration count in
  an iterative solver, so "the same up to rounding" is no guarantee worth a mode; three remain.
- The fastest by default with the strict modes for the tests: the defaults keep today's behaviour (`None` is it), so
  the papers' iteration counts and ADR-0008 stand; the tests set the mode they check.
- MUMPS's determinism options: its FAQ documents determinism only for the parallel analysis (SCOTCH, ParMETIS builds),
  not used here; `ICNTL(28)` changed nothing (measured).
- Two axes (threads, processes): the results already depend on the number of processes by design (ADR-0011), so one
  axis, the threads, with the number of processes fixed.
