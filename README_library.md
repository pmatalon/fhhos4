# fhhos4: the library

The shared library `fhhos4` offers an algebraic multigrid method for the solution of linear systems 
arising from hybrid discretizations (HHO, HDG...) of scalar diffusion and convection-diffusion problems, 
of the form
```
[ A_TT  A_TF ] [x_T]   [b_T]
[ A_FT  A_FF ] [x_F] = [b_F]
```
While the smoothers work on the condensed system (given by the elimination of the cell unknowns)
```
A = A_FF - A_FT A_TT^{-1} A_TF
b = b_F  - A_FT A_TT^{-1} b_T
```
the aggregation method used for the construction of coarse levels requires the blocks of the uncondensed matrix
(to recover the cell-to-cell connectivity amongst other things). Consequently, instead of giving only the condensed system like any black-box AMG method, 
the blocks of the uncondensed matrix must in our case be passed to the solver.

The multigrid hierarchy unfolds as follows:
First, the high polynomial orders are handled by p-multigrid levels, successively reducing the order to the lowest one.
Then, once the system is reduced to one unknown per cell and face, h-coarsening follows, based on the algorithm described in

> D. A. Di Pietro, F. Hülsemann, P. Matalon, P. Mycek, U. Rüde, *Algebraic multigrid preconditioner for statically
> condensed systems arising from lowest-order hybrid discretizations*, SIAM J. Sci. Comput., 2023.

## Applicability

The solver typically applies to symmetric positive definite systems (e.g. diffusion)
and non-symmetric ones (e.g. convection-diffusion, both diffusion- and convection-dominated).
The polynomial bases used for the discretization **must be hierarchical** and ordered by degree for the p-multigrid part to work.

## How to use

Two types of user interface are proposed. 
- (Recommended) If the matrix is symmetric (resp. unsymmetric), pass the three (resp. four) blocks of the uncondensed matrix, namely `A_TT, A_TF, A_FF` (resp. +  `A_FT`)
and the two blocks of the right-hand side `b_T, b_F`.
Through this interface, the condensed system is built internally and the system is solved for both the face and cell unknowns `x_T, x_F`.
- The condensed system can also directly be passed to the solver (along with the blocks). The solver won't perform the condensation itself and solve for the face unknowns only.

The second interface was made solely to make the integration of the `fhoos4` solver less involved in existing codes, which already assemble the condensed matrix and right-hand side, and recover the cell unknowns at the after solving the system.
This way is less efficient than using the first interface, because both the external code and the `fhoos4` solver will invert the `A_TT` block.

| Interface | Setup | Solve | For |
|---|---|---|---|
| **Uncondensed system** | `SetupSymmetric()`, `SetupNonsymmetric()` | `Solve(b_T, b_F, x_T, x_F)` | New codes: assemble the hybrid system, nothing else. The library condenses it and solves it for the face and cell unknowns |
| **Condensed system** | `SetupSymmetricCondensed()`, `SetupNonsymmetricCondensed()` | `Solve(b, x)` | Codes that already condense the matrix: least integration effort |

Symmetric positive definite systems, and non-symmetric ones with a positive definite symmetric part
(convection-diffusion). Krylov method included (FCG, GCR, BiCGSTAB). Parallel (OpenMP).

## Install from source

The library alone needs a C++17 compiler, CMake (>= 3.20), Ninja or Make, and
[Eigen](https://eigen.tuxfamily.org).

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

### The uncondensed system (recommended)

```cpp
#include <fhhos4>

fhhos4::Solver solver;
solver.Dimension = 3;  // to set
solver.FaceDegree = k; // to set
solver.CellDegree = k; // to set
solver.SetupSymmetric(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne);
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
solver.SetupSymmetricCondensed(A, A_TT, A_TF, cellInterpOfOne, faceInterpOfOne);
Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
fhhos4::Result result = solver.Solve(b, x); // the condensed system
```

**The solver keeps a reference to `A`**: `A` must stay alive and unchanged until the next setup or the destruction of
the solver (a temporary `A` does not compile). The blocks may be freed after the setup. A new code should use the
uncondensed system.

### Non-symmetric matrices

For a non-symmetric matrix (e.g. convection-diffusion), give also the block `A_FT` (faces x cells):
`SetupNonsymmetric(A_TT, A_TF, A_FT, A_FF, ...)` then `Solve(b_T, b_F, x_T, x_F)`, or
`SetupNonsymmetricCondensed(A, A_TT, A_TF, A_FT[, A_FF], ...)` then `Solve(b, x)`. The same algorithm runs on the
matrix itself, with GCR instead of FCG and Gauss-Seidel smoothers that follow the direction of the flow (downwind
order).

### Inputs

- **The blocks** `A_TT` (cells x cells), `A_TF` (cells x faces), `A_FF` (faces x faces), and `A_FT` (faces x cells) for
  a non-symmetric matrix. The DoFs of each cell (resp. face) are contiguous, cells and faces in any order. In a
  symmetric system, `A_TT` and `A_FF` are read from their lower triangular part. With the condensed interface, `A`
  (both triangles), and `A_FF` is optional.
- **Matrix types**: `Eigen::SparseMatrix<double>`, column- or row-major, `int` or `long` indices (all of the same
  type). Matrices of type `fhhos4::SparseMatrix` (row-major, fhhos4's index type), compressed, are read in place;
  the others are copied.
- **The interpolation of the function 1** on the cell and face bases: one coefficient per cell (`cellInterpOfOne`, in
  the order of the blocks of rows of `A_TF`) and per face (`faceInterpOfOne`, blocks of columns). The bases must be
  **hierarchical with a constant first function** `phi_0`: `1 = c phi_0`, i.e. `c = 1/phi_0`, which may vary from one
  cell or face to another (e.g. L2-orthonormal bases: `c = sqrt(|T|)` on a cell `T`); `c = 1` everywhere if `phi_0 = 1`.
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
