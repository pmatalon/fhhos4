# fhhos4_AMG: the AMG of fhhos4 for other codes

The shared library `fhhos4_AMG` (short name: fhhos4) gives other codes the algebraic multigrid of fhhos4 for the
statically condensed systems of hybrid discretizations (HHO, HDG...) of scalar diffusion problems. It is the
Uncondensed AMG (U-AMG) of

> D. A. Di Pietro, F. Hülsemann, P. Matalon, P. Mycek, U. Rüde, *Algebraic multigrid preconditioner for statically
> condensed systems arising from lowest-order hybrid discretizations*, SIAM J. Sci. Comput., 2023.

extended to k >= 1: it first coarsens the polynomial degree down to k=0 (p-levels), then aggregates the cells and faces
as in the paper (h-levels).

It solves the hybrid system of cell and face unknowns by static condensation: the AMG solves the condensed system
`A x_F = b` (the face unknowns). Unlike a black-box AMG, it builds its coarse levels from the blocks of the uncondensed
system:

```
[ A_TT  A_TF ] [x_T]   [b_T]   cells        A = A_FF - A_FT A_TT^{-1} A_TF
[ A_FT  A_FF ] [x_F] = [b_F]   faces        b = b_F  - A_FT A_TT^{-1} b_T
```

## Build and link

The library is built with fhhos4 (target `fhhos4_AMG`, `build/lib/libfhhos4_AMG.so`, see the main README). Alone,
it only needs Eigen, CMake and a C++17 compiler (no CGAL, GMSH, GoogleTest; ~45 s and 1.5 GB of RAM):

```bash
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_PROGRAM=OFF ..
```

Its calling code only needs Eigen, the same version as fhhos4's (see below). With CMake:

```cmake
find_package(fhhos4 REQUIRED)
target_link_libraries(my_code PRIVATE fhhos4::fhhos4)
```

and `-Dfhhos4_DIR=<fhhos4>/build` (from the build tree), or, after `cmake --install build --prefix <prefix>
[--component library]`, `-DCMAKE_PREFIX_PATH=<prefix>`. CMake puts the directory of the library in the RPATH of the executable.

## Version

- `find_package(fhhos4 1.2 REQUIRED)` accepts the versions 1.x with x >= 2: the API stays source-compatible within a
  major version.
- Binary compatibility only within a minor version (`libfhhos4_AMG.so.1.2`): adding a parameter to `fhhos4::Solver`
  changes its size. Rebuild the calling code when fhhos4 changes minor version.
- `FHHOS4_VERSION` (`"1.2.0"`), `FHHOS4_VERSION_MAJOR`, `_MINOR`, `_PATCH`: the version of the header (`<fhhos4>`).
- `fhhos4::Version()`: the version of the library, followed by the commit (`git describe`) when it was built from a git
  clone, e.g. `1.3.0 (v1.3.0-12-gabc1234-dirty)`. Give it in bug reports.

## Use

### The uncondensed system: the interface of choice

```cpp
#include <fhhos4> // or <fhhos4.h>

fhhos4::Solver solver;
solver.Dimension = 3;  // to set
solver.FaceDegree = k; // to set
solver.CellDegree = k; // to set
solver.SetupSymmetricFromBlocks(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne);
Eigen::VectorXd x_T(b_T.size()), x_F = Eigen::VectorXd::Zero(b_F.size());
fhhos4::Result result = solver.Solve(b_T, b_F, x_T, x_F); // flexible conjugate gradient preconditioned by the AMG
```

The calling code assembles the uncondensed system, without the DoFs of the Dirichlet faces (moved to the right-hand
side: `b_T -= A_TD x_D`, `b_F -= A_FD x_D`), and nothing else: no static condensation. The library condenses the system,
solves the condensed system and computes the cell unknowns, `x_T = A_TT^{-1} (b_T - A_TF x_F)`. The condensed matrix is
available for other uses: `CondensedMatrix()`. `SetupSymmetricFromBlocks()` (symmetric positive definite system)
builds the levels once: call `Solve()` for each right-hand side.

The library condenses by cells, in parallel: on HArDCore3D's 32^3 system at k = 1 (285,696 face DoFs, 9.2 million
non-zeros in `A`), the setup from the blocks takes 0.3 s more than the setup from `A` on 1 thread (0.71-0.80 s instead of
0.38-0.50 s), 0.08 s more on 8 threads (0.25-0.35 s instead of 0.17-0.29 s); the solve is the same (see
`PERFORMANCE.md`). It computes the rows of `A` in the pattern of `A_FF`, completed (and sorted) only in the rows where
`A_FF` lacks couplings of the faces of a cell (zero coefficients dropped by an assembly): 1.3x (1 thread) to 1.8x (8
threads) slower when all the rows are completed.

### The condensed system: to integrate the library into existing codes

For a code that condenses statically already (e.g. HArDCore3D and Disk++, in their assembly loops), the least effort:
give the library the condensed matrix with the blocks, and keep the code's own condensation of the right-hand side and
computation of the cell unknowns:

```cpp
solver.SetupSymmetric(A, A_TT, A_TF, cellInterpOfOne, faceInterpOfOne);
Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
fhhos4::Result result = solver.Solve(b, x); // the condensed system
```

These setups exist only for that: a new code should use the uncondensed system.

### Non-symmetric matrices

For a **non-symmetric** matrix (e.g. convection-diffusion, whose symmetric part is positive definite), give also the
block `A_FT` (faces x cells): `SetupNonsymmetricFromBlocks(A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne)`
then `Solve(b_T, b_F, x_T, x_F)`, or `SetupNonsymmetric(A, A_TT, A_TF, A_FT[, A_FF], ...)` then `Solve(b, x)`, all the
blocks read in full. The same algorithm runs, on
the matrix itself (no symmetric part is formed); only the optimizations of the symmetric case are off (reading the
lower triangular parts, Cholesky factorizations; the local Galerkin products compute full matrices):
- the couplings of the aggregation are computed from the averages of the traces of the blocks `A_TF(T, F)` and
  `A_FT(F, T)` (equal for a symmetric matrix). On HArDCore3D's convection-diffusion, the couplings of the rows of the
  matrix (as AGMG's strong couplings: from `A_TF(T1, F)` and `A_FT(F, T2)`) cost 20% more iterations in the
  diffusion-dominated case at k = 0, and more without the downwind order;
- the reconstruction of the cells from the faces is `-A_TT^{-1} A_TF`, the inverse of the static condensation (LU of
  the blocks of `A_TT`), and the smoothing of the prolongations uses the rows of the matrix;
- the coarse operators are the Galerkin products `P^T A P` (restriction `P^T`);
- the Gauss-Seidel smoothers of each level follow the **downwind order** of its operator, in the pre- and
  post-smoothing (Bey & Wittum 1997, Hackbusch & Probst 1997): the faces are sorted so that each one comes after those
  it is strongly and one-sidedly coupled to (convection), cycles (recirculation) broken. In an order unrelated to the
  flow, a sweep carries the information only a few cells downstream, and the iterations grow with the number of cells
  along the flow (HArDCore3D, convection-dominated, 32^3, k = 1: 16 iterations, 4 with the downwind order; see
  `PERFORMANCE.md`). Without one-sided coupling (diffusion-dominated), the natural order is kept. On a level with
  one-sided couplings, the configured direction of the Gauss-Seidel smoothers is not respected: all their sweeps,
  backward and symmetric ones included, follow the downwind order (Jacobi and SOR ignore it);
- the K-cycle iterates with GCR, the default Krylov method is GCR and the default coarse solver LU (FCG and Cholesky
  are refused).

### Krylov method

The Krylov method is part of the library (`Krylov`: `"fcg"` by default for a symmetric matrix, `"gcr"`, `"bicgstab"`,
or `"none"` for the cycles alone), which does not give the multigrid cycle alone. With the default K-cycle, the AMG is a
**non-linear preconditioner** (it solves the coarse levels with inner FCG or GCR iterations), so it requires a flexible
Krylov method: the flexible conjugate gradient (FCG), or GCR (generalized conjugate residual, restarted every 10
iterations, in the economical version of Notay, ETNA 37, 2010, Algorithm 3.1; it minimizes the residual norm: for any
matrix; the K-cycle then iterates with GCR too). BiCGSTAB is not flexible: it requires a linear cycle (`Cycle = 'V'` or
`'W'`), and the setup refuses the K-cycle with it.

### Inputs

- **The matrices**, without the DoFs of the Dirichlet faces: the blocks of the uncondensed system `A_TT` (cells x
  cells), `A_TF` (cells x faces), `A_FF` (faces x faces), and `A_FT` (faces x cells) for a non-symmetric matrix, which the
  calling code assembles from its local matrices (the rows of `A_TT` and `A_TF` each belong to one cell). The DoFs of
  each cell (resp. face) are contiguous, cells and faces in any order. In a symmetric system, `A_TT` and `A_FF` are
  read from their lower triangular part. The setups from the blocks keep what `Solve(b_T, b_F, x_T, x_F)` needs (the
  inverses of the blocks of `A_TT`, `A_TT^{-1} A_TF`, `A_FT`): the blocks may be freed after the setup.
  `Eigen::SparseMatrix<double>`, column- or row-major, `int` or `long` indices (all of the same type). Of type
  `fhhos4::SparseMatrix` (row-major, with fhhos4's index type: `int`, `intmax_t` if compiled with `-DSMALL_INDEX=OFF`)
  and compressed, they are read in place, without copy; otherwise, they are converted (copied).
- **The condensed matrix `A`** (faces x faces, both triangles), with the condensed interface only: `SetupSymmetric(A,
  A_TT, A_TF[, A_FF], ...)`; `A_FF` is then optional, as the algorithm does not use it. **The solver keeps a reference to
  `A`**: `A` must stay alive and unchanged until the next setup or the destruction of the solver (a temporary `A` does
  not compile). The blocks are only read by the setup: they may be freed after it. In a symmetric system, the FCG
  applies `A` from its lower triangular part: for an `A` symmetric up to rounding (e.g. assembled from local Schur
  complements), the solution is that of its lower triangle, which differs from a direct solver's on the full `A` by
  ~cond(A) times the asymmetry (Disk++'s 2D systems at k = 2, asymmetry 3e-16: up to 1% of L2 errors of 3e-11).
- **The right-hand side**: `b_T` (one coefficient per row of `A_TF`) and `b_F` (per column), or the condensed `b`.
- **The interpolation of the function 1** on the cell and face bases (`cellInterpOfOne`: one coefficient per row of
  `A_TF`; `faceInterpOfOne`: per column). The bases must be **hierarchical with a constant first function** `phi_0`:
  the interpolation of 1 is then `1/phi_0` on the first DoF of each cell or face, 0 on the others. It may vary (e.g.
  `sqrt|T|` for L2-orthonormal bases, whose constant is `1/sqrt|T|`): the AMG rescales the constant modes before
  coarsening. Without this rescaling, L2-orthonormal bases slow the AMG down (Cube-tet n=16, k=0: 44 iterations instead
  of 26) or make it fail (graded meshes).
- **The dimension and the degrees**: `Dimension`, `FaceDegree` (k) and `CellDegree` (l) must be set. They give the
  DoFs per face and per cell (full polynomial spaces: `C(k+d-1, d-1)` and `C(l+d, d)`), against which the setup checks
  the sizes and structure of the inputs (block-diagonal `A_TT`, interpolation of 1 non-zero only on the first DoF of
  each cell and face).

### Parameters

The public members of `fhhos4::Solver` (see `include/fhhos4.h`), read by the setup (`Tolerance` and `MaxIterations` by
`Solve()`). The defaults are those of the program fhhos4 (`-s fcguamg`), used in the papers: FCG with tolerance 1e-8,
K-cycle, one block Gauss-Seidel iteration before and after, coarsening factor 3.8, Cholesky on the coarsest level
(<= 1000 rows). Except the smoothers on several threads: by default (`PreSmoother` and `PostSmoother` empty),
the setup then chooses the hybrid block Gauss-Seidel (`hbgs`, `hrbgs`: Gauss-Seidel in the rows of each thread, Jacobi
between the threads: a row reads the values of the other threads' rows from the start of the sweep; with the scalar
kernels at k = 0), parallel, unlike the block Gauss-Seidel. Its results depend on the number of threads, not on the
run. On a level where the rows of the threads are poorly separated (more than 20% of the non-zeros couple rows of
different threads: a numbering without locality, such as GMSH's on fhhos4's tetrahedral meshes), it runs the
sequential sweeps, as on 1 thread. It takes a few more iterations than the block Gauss-Seidel (0-6 with FCG, a few
more with BiCGSTAB), but solves faster on several threads: with FCG on 8 threads, 15-20% at k = 0, 35-45% at k = 1, 2
(HArDCore3D's cubes and Voronoi mesh); with BiCGSTAB or the multigrid alone, 2-45% (measured with an earlier version,
whose iterations varied from run to run; see `PERFORMANCE.md`). With GCR (non-symmetric matrices), HArDCore3D's
diffusion-dominated convection-diffusion at 48^3, k = 1: 27 iterations on 4 threads, 26 on 1, solve 3.3-3.5 s instead
of 5.8 s.
The FCG also replaces its residual (updated by the recurrence) by the true residual `b - A x` every 20 iterations and
when it meets the tolerance (GCR: at each restart and when it meets the tolerance), then continues if the true
residual does not: the tolerance applies to the true residual
(`Result::RelativeResidual`). Near the attainable accuracy (~1e-11 relative on Disk++'s 2D systems at k = 2), the
recurrence drifts from the true residual (8e-13 for a true 1.4e-11). A tolerance below the attainable accuracy is
never met: `Solve()` then runs `MaxIterations` iterations, not converged.
`Threads` sets the number of OpenMP threads during the calls (default: the current OpenMP setting, but at most one
thread per physical core unless `OMP_NUM_THREADS` or OpenMP's binding is set: the hyperthreads slow the AMG down).
`Verbosity`: 0 (default) prints nothing.

## Constraints

- **Symmetric positive definite** condensed system (`SetupSymmetric...()`), or non-symmetric with a positive definite
  symmetric part (`SetupNonsymmetric...()`): not the pure Neumann problem (singular).
- **Same Eigen**: Eigen objects cross the boundary of the library, so the calling code must use the same version of
  Eigen and the same alignment settings (the same `-march` family: with or without AVX). The constructor of
  `fhhos4::Solver` checks it and throws otherwise. fhhos4 is built with the conda compilers' `-march=nocona`: compile the calling code
  without `-march=native`, or rebuild fhhos4 with `-DENABLE_NATIVE_ARCH=ON`.
- **Errors**: `fhhos4::Error` (a `std::runtime_error`). `Solve()` does not throw if the tolerance is not reached: see
  `Result::Converged`.
- **Not thread-safe**: one `fhhos4::Solver` object must not be used by several threads at the same time. The library
  parallelizes its own work with OpenMP (on 1 thread, its solve phase is sequential Gauss-Seidel smoothing).

## Tests

`tests/UAMGLibraryTest.cpp`: the program runs U-AMG through the library (`-s libuamg`, `fcglibuamg`, `gcrlibuamg`),
through both interfaces (the uncondensed system by default, the condensed one with `-libuamg-condensed`), and must get
the results of the U-AMG compiled in it (`uamg`, `fcguamg`, `gcruamg`); the public API alone on a hybrid system built in
the test (the uncondensed system against a direct solver, matrices computed by the library, scaled bases, GCR,
BiCGSTAB, a non-symmetric system, checks of the parameters and inputs, silence at verbosity 0); the non-symmetric setup
on the program's symmetric matrices (`-uamg-nonsym`) gives the iterations of the symmetric one.
`tests/StaticCondensationTest.cpp`: the condensation by cells gives the matrices of the global products, bit for bit.

## License

LGPL-3.0-or-later, like the rest of fhhos4 (see the main README). The library only depends on Eigen (MPL-2.0), not on
the GPL dependencies of the program (GMSH, CGAL): codes under other licenses can link it, under the conditions of
section 4 of the LGPL (in short: keep the license notices, and let the user replace the shared library).
