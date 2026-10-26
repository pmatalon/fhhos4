# Block Chebyshev, the library's default smoother for symmetric systems

For a symmetric matrix the library smooths with block Chebyshev of the fourth kind, degree 2, on any number of threads
(2026-10-15): its results do not depend on the number of threads, and it is faster than the hybrid block Gauss-Seidel
on 8 threads in most cases (up to 2x on GMSH tetrahedra at k >= 1), slower on 1 thread. The program keeps the papers'
block Gauss-Seidel (`bgs,rbgs`; ADR-0008).
