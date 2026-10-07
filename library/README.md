# fhhos4_AMG: the AMG of fhhos4 for other codes

The shared library `fhhos4_AMG` (short name: fhhos4) gives other codes the algebraic multigrid of fhhos4 for the
statically condensed systems of hybrid discretizations (HHO, HDG...) of scalar diffusion problems. It is the
Uncondensed AMG (U-AMG) of

> D. A. Di Pietro, F. Hülsemann, P. Matalon, P. Mycek, U. Rüde, *Algebraic multigrid preconditioner for statically
> condensed systems arising from lowest-order hybrid discretizations*, SIAM J. Sci. Comput., 2023.

extended to k >= 1: it first coarsens the polynomial degree down to k=0 (p-levels), then aggregates the cells and faces
as in the paper (h-levels).

It solves the condensed system `A x = b` (the face unknowns). Unlike a black-box AMG, it builds its coarse levels from
the blocks of the uncondensed system, before static condensation:

```
[ A_TT  A_TF ] [x_T]   cells
[ A_FT  A_FF ] [x_F]   faces          A = A_FF - A_FT A_TT^{-1} A_TF
```

## Build and link

The library is built with fhhos4 (target `fhhos4_AMG`, `build/lib/libfhhos4_AMG.so`, see the main README). Its
calling code only needs Eigen, the same version as fhhos4's (see below). With CMake:

```cmake
find_package(fhhos4 REQUIRED)
target_link_libraries(my_code PRIVATE fhhos4::AMG)
```

and `-Dfhhos4_DIR=<fhhos4>/build` (from the build tree), or, after `cmake --install build --prefix <prefix>`,
`-DCMAKE_PREFIX_PATH=<prefix>`. CMake puts the directory of the library in the RPATH of the executable.

## Use

```cpp
#include <fhhos4>

fhhos4::Solver solver;
solver.Dimension = 3;  // to set
solver.FaceDegree = k; // to set
solver.Setup(A, A_TT, A_TF, cellInterpOfOne, faceInterpOfOne);
Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
fhhos4::Result result = solver.Solve(b, x); // flexible conjugate gradient preconditioned by the AMG
```

`Setup()` builds the levels once: call `Solve()` for each right-hand side.

The Krylov method is part of the library (`Krylov`: `"fcg"` by default, `"bicgstab"`, or `"none"` for the cycles
alone), which does not give the multigrid cycle alone. With the default K-cycle, the AMG is a **non-linear
preconditioner** (it solves the coarse levels with inner FCG iterations), so it requires a flexible Krylov method: the
flexible conjugate gradient. BiCGSTAB is not flexible: it requires a linear cycle (`Cycle = 'V'` or `'W'`), and the
setup refuses the K-cycle with it.

### Inputs

- **The matrices**, without the DoFs of the Dirichlet faces: the condensed matrix `A` (faces x faces, both triangles)
  and the blocks of the uncondensed system `A_TT` (cells x cells) and `A_TF` (cells x faces), which the calling code
  assembles from its local matrices (the rows of `A_TT` and `A_TF` each belong to one cell). `A_FF` (faces x faces) is
  optional, as the algorithm does not use it: `Setup(A, A_TT, A_TF, A_FF, ...)`. Without `A`,
  `SetupFromBlocks(A_TT, A_TF, A_FF, ...)` computes it. The DoFs of each cell (resp. face) are contiguous, cells and
  faces in any order. `A_TT` and `A_FF` are read from their lower triangular part. `Eigen::SparseMatrix<double>`,
  column- or row-major, `int` or `long` indices (all of the same type); they are copied, and may be freed after
  `Setup()`.
- **The interpolation of the function 1** on the cell and face bases (`cellInterpOfOne`: one coefficient per row of
  `A_TF`; `faceInterpOfOne`: per column). The bases must be **hierarchical with a constant first function** `phi_0`:
  the interpolation of 1 is then `1/phi_0` on the first DoF of each cell or face, 0 on the others. It may vary (e.g.
  `sqrt|T|` for L2-orthonormal bases, whose constant is `1/sqrt|T|`): the AMG rescales the constant modes before
  coarsening. Without this rescaling, L2-orthonormal bases slow the AMG down (Cube-tet n=16, k=0: 44 iterations instead
  of 26) or make it fail (graded meshes).
- **The dimension and the degrees**: `Dimension` and `FaceDegree` (k) must be set, `CellDegree` defaults to k. They
  give the DoFs per face and per cell (full polynomial spaces: `C(k+d-1, d-1)` and `C(l+d, d)`), against which
  `Setup()` checks the sizes and structure of the inputs (block-diagonal `A_TT`, interpolation of 1 non-zero only on the
  first DoF of each cell and face).

### Parameters

The public members of `fhhos4::Solver` (see `include/fhhos4`), read by `Setup()` (`Tolerance` and `MaxIterations` by
`Solve()`). The defaults are those of the program fhhos4 (`-s fcguamg`), used in the papers: FCG with tolerance 1e-8,
K-cycle, one block Gauss-Seidel iteration before and after, coarsening factor 3.8, Cholesky on the coarsest level
(<= 1000 rows). `Threads` sets the number of OpenMP threads during the calls (default: the current OpenMP setting).
`Verbosity`: 0 (default) prints nothing.

## Constraints

- **Symmetric positive definite** condensed system: not the pure Neumann problem (singular).
- **Same Eigen**: Eigen objects cross the boundary of the library, so the calling code must use the same version of
  Eigen and the same alignment settings (the same `-march` family: with or without AVX). The constructor of
  `fhhos4::Solver` checks it and throws otherwise. fhhos4 is built with the conda compilers' `-march=nocona`: compile the calling code
  without `-march=native`, or rebuild fhhos4 with `-DENABLE_NATIVE_ARCH=ON`.
- **Errors**: `fhhos4::Error` (a `std::runtime_error`). `Solve()` does not throw if the tolerance is not reached: see
  `Result::Converged`.
- **Not thread-safe**: one `fhhos4::Solver` object must not be used by several threads at the same time. The library
  parallelizes its own work with OpenMP; its solve phase is mostly sequential (Gauss-Seidel smoothing).

## Tests

`tests/UAMGLibraryTest.cpp`: the program runs U-AMG through the library (`-s libuamg`, `fcglibuamg`) and must get the
results of the U-AMG compiled in it (`uamg`, `fcguamg`); the public API alone on a hybrid system built in the test
(matrices computed by the library, scaled bases, BiCGSTAB, checks of the parameters and inputs, silence at verbosity 0).
