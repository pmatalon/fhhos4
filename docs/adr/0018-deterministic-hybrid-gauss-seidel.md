# Deterministic hybrid Gauss-Seidel, sequential where threads share too much

For a non-symmetric matrix on several threads the default smoothers are hybrid block Gauss-Seidel (`hbgs`, `hrbgs`):
Gauss-Seidel in the rows of each thread, Jacobi between the threads through a copy of the values the others read, so
the results depend on the number of threads but not on the run (2026-10-12). A level where more than 20% of the
non-zeros couple rows of different threads sweeps sequentially.

## Considered Options

- The racy version (threads reading each other's values as they change): with GCR, 191 iterations instead of 26 on
  8 threads (HArDCore3D 48^3, k = 1).
- Jacobi between the threads everywhere: diverged on GMSH's tetrahedra (32% of the non-zeros between threads).
- l1-Gauss-Seidel (Baker et al. 2011): tetrahedra 35/72 iterations instead of 28/33, convection 5 -> 10.
