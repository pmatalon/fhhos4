# Deterministic hybrid Gauss-Seidel, in l1 Gauss-Seidel by default, sequential where threads share too much

For a non-symmetric matrix on several threads the smoothers are hybrid block Gauss-Seidel: Gauss-Seidel in the rows of
each thread, Jacobi between the threads through a copy of the values the others read, so the results depend on the
number of threads but not on the run (2026-10-12). A level where more than 20% of the non-zeros couple rows of
different threads sweeps sequentially.

Since 2026-11-02 the default is the l1 hybrid block Gauss-Seidel (`L1HybridBlockGaussSeidel`, `l1hbgs`), a smoother of
its own: the rows coupled to other threads and processes in l1 Gauss-Seidel (theta = 1/2: half the sum of the |a_ij| of
these couplings on their diagonal), convergent whatever the parts. Under the 20% threshold, the plain hybrid sweeps
diverged at k = 2 in 3D (Cube-tet n=8 k=2, 2 threads: 200 iterations, 34 on one; l1: 94), and so did Jacobi between
processes. Since 2026-11-03 the scalar sweeps of degree 0 too (`GaussSeidel` with `l1`): Jacobi between processes
stagnated there on the finest level (Cube-tet n=32 k=0, `-uamg-nonsym`, 2 and 4 processes: GCR at 3.7e-2 after 200
iterations; l1: 29-30, 28 on one process); on the k=0 levels below a p-coarsening it changed nothing on several
processes (k=1 and k=2: the same iterations with and without). Its cost on one process, the l1 terms of the threads'
chunks (HArDCore3D, 4-8 threads): +0-2 iterations, 4 -> 8 at k=0 convection-dominated, times within the noise. The plain hybrid one (`HybridBlockGaussSeidel`) stays, with its former meaning: where it converges it takes
fewer iterations (HArDCore3D, 4 threads: convection-dominated 4 -> 7, 8 -> 11, tetrahedra 9 -> 14 with l1; diffusion
+0-2). A smoother rather than a parameter: an explicit choice keeps its meaning, only the default changes; as hypre
and MFEM, which have l1-Gauss-Seidel as its own relaxation type.

## Considered Options

- The racy version (threads reading each other's values as they change): with GCR, 191 iterations instead of 26 on
  8 threads (HArDCore3D 48^3, k = 1).
- Jacobi between the threads everywhere: diverged on GMSH's tetrahedra (32% of the non-zeros between threads).
- l1-Gauss-Seidel with the full l1 term (theta = 1, Baker et al. 2011, hypre): tetrahedra 35/72 iterations instead of
  28/33, convection 5 -> 10.
- The plain hybrid one by default: diverges at k = 2 in 3D, as above.
- A boolean parameter for l1 (as Trilinos/Ifpack2's `relaxation: use l1`): would change the meaning of an explicit
  `HybridBlockGaussSeidel`.
- The l1 terms only where the diagonal is not large enough (Ifpack2's eta, hypre's truncated variant): not tried yet;
  could recover the iterations where the plain one converges.
- The rows between threads by colours (a true Gauss-Seidel, a barrier per colour): not tried.
