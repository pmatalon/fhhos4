# Deterministic hybrid Gauss-Seidel with an l1 shift, sequential where threads share too much

For a non-symmetric matrix on several threads the default smoothers are hybrid block Gauss-Seidel (`hbgs`, `hrbgs`):
Gauss-Seidel in the rows of each thread, Jacobi between the threads through a copy of the values the others read, so
the results depend on the number of threads but not on the run (2026-10-12). A level where more than 20% of the
non-zeros couple rows of different threads sweeps sequentially.

Since 2026-11-02 the rows coupled to other threads (and processes) have an l1 shift, theta = 1/2 (half the sum of the
|a_ij| of these couplings on their diagonal): under the 20% threshold, Jacobi between the threads diverged at k = 2 in
3D (Cube-tet n=8 k=2, 2 threads: 200 iterations, 34 on one; with the shift 94). It costs where the plain hybrid sweeps
converged: HArDCore3D 32^3 k=1 convection-dominated 5 -> 7 iterations, diffusion-dominated +0-1.

## Considered Options

- The racy version (threads reading each other's values as they change): with GCR, 191 iterations instead of 26 on
  8 threads (HArDCore3D 48^3, k = 1).
- Jacobi between the threads everywhere: diverged on GMSH's tetrahedra (32% of the non-zeros between threads).
- l1-Gauss-Seidel (Baker et al. 2011), the full shift (theta = 1): tetrahedra 35/72 iterations instead of 28/33,
  convection 5 -> 10.
- No shift (until 2026-11-02): diverges at k = 2 in 3D, as above.
- The rows between threads by colours (a true Gauss-Seidel, a barrier per colour): not tried; the next option if the
  shift's iterations cost too much.
