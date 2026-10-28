# The library's defaults differ from the program's

The program's default options stay those of the papers, whose iteration counts the tests check: a new numerical method
(prolongation, precision) is an opt-in option there. One exception: the smoother of U-AMG, block Chebyshev in both
since 2026-10-28 (ADR-0017), the tests setting the papers' Gauss-Seidel. The library takes the fastest validated
methods as its defaults: the block Chebyshev smoother, the harmonic extension as coarsening prolongation (2026-10-22),
mixed precision (ADR-0010), the true residual in the stopping test, and GCR with LU for a non-symmetric system. So the
library's documentation reads its defaults in its own code (`library/include/fhhos4.h`, `library/Solver.cpp`), never
from the program.
