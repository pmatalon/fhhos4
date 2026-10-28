# Block Chebyshev, the default smoother of U-AMG for symmetric systems

For a symmetric matrix U-AMG smooths with block Chebyshev of the fourth kind, degree 2: its results do not depend on the
numbers of threads and processes. In the library since 2026-10-15: faster than the hybrid block Gauss-Seidel on 8
threads in most cases (up to 2x on GMSH tetrahedra at k >= 1), slower on 1 thread. In the program since 2026-10-28,
even on one process, so that a run gives the same results on any number of processes: on several, Gauss-Seidel is
Jacobi between them, which diverged (3D, k=2: FCG stagnated at 200 iterations, 32 on one process; Chebyshev 50 on 1, 2
and 4). The papers' block Gauss-Seidel stays an option (`-smoothers bgs`), which the tests and `reproducibility/` set.

## Considered Options

- l1-Gauss-Seidel between the processes (the diagonal of the rows coupled to other processes shifted by half the sum of
  their off-process |a_ij|: a contraction in the energy norm for any partition): 61 iterations on 2 and 4 processes,
  32 on one; it still depends on the partition.
- Chebyshev on several processes only: the results of the program would change with the number of processes.
