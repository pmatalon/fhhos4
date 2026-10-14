# fhhos4: the library

The shared library `fhhos4` (`libfhhos4.so`) gives your code the algebraic multigrid of fhhos4 for the statically
condensed systems of hybrid discretizations (HHO, HDG...) of scalar diffusion problems: the Uncondensed AMG (U-AMG) of

> D. A. Di Pietro, F. Hülsemann, P. Matalon, P. Mycek, U. Rüde, *Algebraic multigrid preconditioner for statically
> condensed systems arising from lowest-order hybrid discretizations*, SIAM J. Sci. Comput., 2023.

extended to any degree. The AMG solves the condensed system `A x_F = b` (the face unknowns), but, unlike a black-box
AMG, it builds its coarse levels from the blocks of the uncondensed system:

```
[ A_TT  A_TF ] [x_T]   [b_T]   cells        A = A_FF - A_FT A_TT^{-1} A_TF
[ A_FT  A_FF ] [x_F] = [b_F]   faces        b = b_F  - A_FT A_TT^{-1} b_T
```

- **p-levels** (degree k >= 1): the first levels lower the polynomial degree of the faces and cells (by 2 at each
  level) down to 0, by keeping the low-degree modes of the hierarchical bases.
- **h-levels** (degree 0): the cells are aggregated pairwise, each with its most strongly coupled neighbour, the
  strength of the coupling of two neighbouring cells coming from the blocks `A_TF(T, F)` of both cells on their common
  face (several passes per level, until the cells are `CoarseningFactor` times fewer). The faces inside an aggregate
  are removed, those between two aggregates merged into one coarse face. The prolongation reconstructs the cell
  unknowns of each aggregate from its coarse faces (`-A_TT^{-1} A_TF`, the inverse of the static condensation) and
  takes their trace on the removed faces; the coarse operator is the Galerkin product.

## Interfaces

| Interface | Setup | Solve | For |
|---|---|---|---|
| **Uncondensed system** | `SetupSymmetricFromBlocks()`, `SetupNonsymmetricFromBlocks()` | `Solve(b_T, b_F, x_T, x_F)` | New codes: assemble the hybrid system, nothing else. The library condenses it (in parallel) and computes the cell unknowns |
| **Condensed system** | `SetupSymmetric()`, `SetupNonsymmetric()` | `Solve(b, x)` | Codes that condense statically already: least integration effort |

Symmetric positive definite systems, and non-symmetric ones with a positive definite symmetric part
(convection-diffusion). Krylov method included (FCG, GCR, BiCGSTAB). Parallel (OpenMP).

## Install from source

The library alone needs a C++17 compiler, CMake (>= 3.20), Ninja or Make, and
[Eigen](https://eigen.tuxfamily.org) (no CGAL, GMSH, GoogleTest; ~45 s and 1.5 GB of RAM). With conda,
`conda env create -f conda/environment.yml` provides all of them; otherwise install Eigen with your package manager.

```bash
git clone https://github.com/pmatalon/fhhos4.git && cd fhhos4
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_PROGRAM=OFF
cmake --build build
cmake --install build --prefix <prefix>    # lib/, include/ and the CMake package
```

Without installing, the library is `build/lib/libfhhos4.so`.

Build options (`-D<OPTION>=ON|OFF`): `ENABLE_OPENMP` (ON: parallel), `ENABLE_NATIVE_ARCH` (OFF: `-march=native`, see the
constraints below), `SMALL_INDEX` (ON: `int` indices in `fhhos4::SparseMatrix`, `intmax_t` otherwise).

## Link

Your code needs the same Eigen as the library (see the constraints below). With CMake:

```cmake
find_package(fhhos4 REQUIRED)
target_link_libraries(my_code PRIVATE fhhos4::fhhos4)
```

and `-DCMAKE_PREFIX_PATH=<prefix>` (installation), or `-Dfhhos4_DIR=<fhhos4>/build` (build tree). CMake puts the
directory of the library in the RPATH of the executable.

## Use

The public interface is `library/include/fhhos4.h` (`#include <fhhos4>`), documented in its comments.

### The uncondensed system: the interface of choice

```cpp
#include <fhhos4>

fhhos4::Solver solver;
solver.Dimension = 3;  // to set
solver.FaceDegree = k; // to set
solver.CellDegree = k; // to set
solver.SetupSymmetricFromBlocks(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne);
Eigen::VectorXd x_T(b_T.size()), x_F = Eigen::VectorXd::Zero(b_F.size());
fhhos4::Result result = solver.Solve(b_T, b_F, x_T, x_F); // flexible conjugate gradient preconditioned by the AMG
```

Your code assembles the uncondensed system, without the DoFs of the Dirichlet faces (moved to the right-hand side:
`b_T -= A_TD x_D`, `b_F -= A_FD x_D`), and nothing else: no static condensation. The library condenses the system,
solves the condensed system and computes the cell unknowns, `x_T = A_TT^{-1} (b_T - A_TF x_F)`. The condensed matrix is
available for other uses: `CondensedMatrix()`. The setup builds the levels once: call `Solve()` for each right-hand
side. The blocks may be freed after the setup.

### The condensed system: to integrate the library into existing codes

For a code that condenses statically already (in its assembly loops), the least effort: give the library the condensed
matrix with the blocks, and keep the code's own condensation of the right-hand side and computation of the cell
unknowns:

```cpp
solver.SetupSymmetric(A, A_TT, A_TF, cellInterpOfOne, faceInterpOfOne);
Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
fhhos4::Result result = solver.Solve(b, x); // the condensed system
```

**The solver keeps a reference to `A`**: `A` must stay alive and unchanged until the next setup or the destruction of
the solver (a temporary `A` does not compile). The blocks may be freed after the setup. A new code should use the
uncondensed system.

### Non-symmetric matrices

For a non-symmetric matrix (e.g. convection-diffusion), give also the block `A_FT` (faces x cells):
`SetupNonsymmetricFromBlocks(A_TT, A_TF, A_FT, A_FF, ...)` then `Solve(b_T, b_F, x_T, x_F)`, or
`SetupNonsymmetric(A, A_TT, A_TF, A_FT[, A_FF], ...)` then `Solve(b, x)`. The same algorithm runs on the matrix itself,
with GCR instead of FCG and Gauss-Seidel smoothers that follow the direction of the flow (downwind order).

### Inputs

- **The blocks** `A_TT` (cells x cells), `A_TF` (cells x faces), `A_FF` (faces x faces), and `A_FT` (faces x cells) for
  a non-symmetric matrix. The DoFs of each cell (resp. face) are contiguous, cells and faces in any order. In a
  symmetric system, `A_TT` and `A_FF` are read from their lower triangular part. With the condensed interface, `A`
  (both triangles), and `A_FF` is optional.
- **Matrix types**: `Eigen::SparseMatrix<double>`, column- or row-major, `int` or `long` indices (all of the same
  type). Matrices of type `fhhos4::SparseMatrix` (row-major, fhhos4's index type), compressed, are read in place;
  the others are copied.
- **The interpolation of the function 1** on the cell and face bases (`cellInterpOfOne`: one coefficient per row of
  `A_TF`; `faceInterpOfOne`: per column). The bases must be **hierarchical with a constant first function** `phi_0`:
  the interpolation of 1 is then `1/phi_0` on the first DoF of each cell or face, 0 on the others. `phi_0` may vary
  from one cell or face to another (e.g. L2-orthonormal bases).
- **The dimension and the degrees**: `Dimension` (2 or 3), `FaceDegree` (k) and `CellDegree` (l) must be set. The setup checks
  the sizes and structure of the inputs against them (full polynomial spaces).

### Parameters

The public members of `fhhos4::Solver`, read by the setup (`Tolerance` and `MaxIterations` by `Solve()`), with the
library's defaults (`library/include/fhhos4.h`, which details each one):

| Parameter | Default | Values |
|---|---|---|
| `Krylov` | `""`: `"fcg"` for a symmetric matrix, `"gcr"` otherwise | `"fcg"` (flexible conjugate gradient), `"gcr"` (GCR(10), any matrix), `"bicgstab"` (V- or W-cycle only), `"none"` (the multigrid cycles alone) |
| `Tolerance` | `1e-8` | `Solve()` stops when the relative residual of the condensed system (the true one, `Result::RelativeResidual`) is below it |
| `MaxIterations` | `200` | |
| `Cycle` | `'K'` | `'V'`, `'W'` or `'K'` (K-cycle: a non-linear preconditioner, which requires FCG or GCR) |
| `PreSmoothingIterations`, `PostSmoothingIterations` | `1`, `1` | |
| `PreSmoother` | `""`: `"hbgs"` on several threads, `"bgs"` on 1 | Smoother codes of the program (`fhhos4 -h`, `-smoothers`). `"bgs"`: block Gauss-Seidel (one block per face); `"hbgs"`: its parallel (hybrid) version |
| `PostSmoother` | `""`: the backward version of `PreSmoother` (`"rbgs"`, `"hrbgs"`) | |
| `CoarseningFactor` | `3.8` | Ratio of the numbers of cells of two successive h-levels |
| `CoarseMatrixMaxSize` | `1000` | The coarsening stops when the matrix has at most this many rows |
| `CoarseSolver` | `""`: `"ch"` (sparse Cholesky) for a symmetric matrix, `"lu"` otherwise | |
| `Threads` | `0`: the OpenMP setting, but at most one thread per physical core unless `OMP_NUM_THREADS` or OpenMP's binding is set | |
| `Verbosity` | `0`: silent | `1`: one line per setup and per solve; `2`: all the messages, with the residual at each iteration |

## Constraints

- **Symmetric positive definite** condensed system (`SetupSymmetric...()`), or non-symmetric with a positive definite
  symmetric part (`SetupNonsymmetric...()`): not the pure Neumann problem (singular).
- **Same Eigen**: Eigen objects cross the boundary of the library, so your code must use the same version of Eigen and
  the same alignment settings (the same `-march` family: with or without AVX). The constructor of `fhhos4::Solver`
  checks it and throws otherwise. The library is built without `-march=native` by default: compile your code without
  it too, or build the library with `-DENABLE_NATIVE_ARCH=ON`.
- **Errors**: `fhhos4::Error` (a `std::runtime_error`). `Solve()` does not throw if the tolerance is not reached: see
  `Result::Converged`.
- **Not thread-safe**: one `fhhos4::Solver` object must not be used by several threads at the same time. The library
  parallelizes its own work with OpenMP.

## License

LGPL-3.0-or-later, like the rest of fhhos4 (see [README.md](README.md)). The library only depends on Eigen (MPL-2.0),
not on the GPL dependencies of the program (GMSH, CGAL): codes under other licenses can link it, under the conditions
of section 4 of the LGPL (in short: keep the license notices, and let the user replace the shared library).

## Version

- `find_package(fhhos4 1.3 REQUIRED)` accepts the versions 1.x with x >= 3: the API stays source-compatible within a
  major version.
- Binary compatibility only within a minor version (`libfhhos4.so.1.3`): adding a parameter to `fhhos4::Solver`
  changes its size. Rebuild your code when fhhos4 changes minor version.
- `FHHOS4_VERSION` (`"1.3.0"`), `FHHOS4_VERSION_MAJOR`, `_MINOR`, `_PATCH`: the version of the header (`<fhhos4>`).
- `fhhos4::Version()`: the version of the library, followed by the commit (`git describe`) when it was built from a git
  clone, e.g. `1.3.0 (v1.3.0-12-gabc1234-dirty)`. Give it in bug reports.
