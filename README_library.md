# the fhhos4 solver as an external library

# 1. Installation

## Install from source

The library alone needs a C++17 compiler, CMake (>= 3.20), Ninja or Make, and
[Eigen](https://eigen.tuxfamily.org) 3.4 or 5.x.

```bash
git clone https://github.com/pmatalon/fhhos4.git && cd fhhos4
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_PROGRAM=OFF
cmake --build build
cmake --install build --prefix <prefix>    # lib/, include/ and the CMake package
```

Without installing, the library is `build/lib/libfhhos4.so`.

Build options (`-D<OPTION>=ON|OFF`): `ENABLE_OPENMP` (ON: parallel), `ENABLE_NATIVE_ARCH` (OFF: `-march=native`:
optimized for the CPU of the build machine, may not run on others), `SMALL_INDEX` (ON: `int` indices in
`fhhos4::SparseMatrix`, `intmax_t` otherwise), `ENABLE_PYTHON` (OFF: the Python module, see [Python](#python)),
`ENABLE_MPI` (OFF: the solver on several processes, see [Several processes](#several-processes-mpi); needs MPI, and
your code is then compiled and linked with MPI too; not with `ENABLE_PYTHON` yet).

### Conda package

The recipe `conda/recipe` builds the conda package `fhhos4` (Linux), with [rattler-build](https://rattler.build)
(`conda install -c conda-forge rattler-build`):

```bash
rattler-build build -r conda/recipe -c conda-forge --experimental   # the package in ./output
conda install -c ./output -c conda-forge fhhos4
```

The package is the library alone. It depends on Eigen only (>= 3.4, see the constraints below), and a package of your
code built with it requires the same minor version of `fhhos4` (see the version below). `find_package(fhhos4)` finds it
in the environment. `--experimental`: the recipe reads the version from `CMakeLists.txt`.

## Link

Your code needs Eigen 3.4 or 5.x (see the constraints below). With CMake:

```cmake
find_package(fhhos4 REQUIRED)
target_link_libraries(my_code PRIVATE fhhos4::fhhos4)
```

and `-DCMAKE_PREFIX_PATH=<prefix>` (installation), or `-Dfhhos4_DIR=<fhhos4>/build` (build tree). CMake puts the
directory of the library in the RPATH of the executable. A library compiled with MPI brings MPI along
(`MPI::MPI_CXX`).

# 2. The algorithm in a nutshell

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
Bases not ordered by degree, but whose first function is constant, work with the p-coarsening straight to the lowest
order (`PCoarsening::ToZero`), or better, reordered by degree: e.g. in the DoF lists of `ExtractSymmetricBlocks()`
(see [Other formats](#other-formats-the-helpers)), at no extra cost.

# 3. How to use

Two types of user interface are proposed. 
- **Interface 1** (Recommended): If the matrix is symmetric (resp. unsymmetric), pass the three (resp. four) blocks of the uncondensed matrix, namely `A_TT, A_TF, A_FF` (resp. +  `A_FT`)
and the two blocks of the right-hand side `b_T, b_F`.
Through this interface, the condensed system is built internally and the system is solved for both the face and cell unknowns `x_T, x_F`.
- **Interface 2**: The condensed system can also directly be passed to the solver (along with the blocks). The solver won't perform the condensation itself and solve for the face unknowns only.

The second interface was made solely to make the integration of the `fhhos4` solver less involved in existing codes, which already assemble the condensed matrix and right-hand side, and recover the cell unknowns after solving the system.
This way is less efficient than using the first interface, because, at the very least, both the external code and the `fhhos4` solver will factorize the `A_TT` block.

For existing codes, the recommendation is to begin with Interface 2 to try out the solver with the least integration effort, and then move on to Interface 1 for better performance. New codes should directly use Interface 1. Interface 1 relieves your code from all matters linked to static condensation, including the recovering of the cell unknowns at the end. 

| Interface | Setup | Solve | For |
|---|---|---|---|
| **Uncondensed system** | `SetupSymmetric()`, `SetupNonsymmetric()` | `Solve(b_T, b_F, x_T, x_F)` | New codes: assemble the hybrid system, nothing else. The library condenses it and solves it for the face and cell unknowns |
| **Condensed system** | `SetupSymmetricCondensed()`, `SetupNonsymmetricCondensed()` | `Solve(b, x)` | Codes that already condense the matrix: least integration effort |

Symmetric positive definite systems, and non-symmetric ones with a positive definite symmetric part
(convection-diffusion). Krylov method included (FCG, GCR, BiCGSTAB). Parallel (OpenMP), and on several processes (MPI,
see [Several processes](#several-processes-mpi)).

### Interface 1 (recommended)

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
`b_T -= A_TD x_D`, `b_F -= A_FD x_D`), and nothing else: no static condensation. The blocks are read in place, in the
CSR format of the library: `fhhos4::SparseMatrix` (Eigen, row-major, compressed), or a `fhhos4::CSRView` of your own
arrays (see [Inputs](#inputs); other formats: the helpers, which copy them). The library condenses the system,
solves the condensed system and computes the cell unknowns, `x_T = A_TT^{-1} (b_T - A_TF x_F)`. The condensed matrix is
available for other uses: `CondensedMatrix()`. The setup builds the levels once: call `Solve()` for each right-hand
side. The blocks may be freed after the setup.

### Interface 2 (less effort for existing codes)

For a code that condenses statically already (in its assembly loops), the least effort: give the library the condensed
matrix with the blocks, and keep the code's own condensation of the right-hand side and computation of the cell
unknowns:

```cpp
solver.SetupSymmetricCondensed(A, A_TT, A_TF, cellInterpOfOne, faceInterpOfOne);
Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
fhhos4::Result result = solver.Solve(b, x); // the condensed system
```

**The solver keeps a reference to `A`** (read in place): `A` must stay alive and unchanged until the next setup or the
destruction of the solver (a temporary `A` does not compile). The blocks may be freed after the setup. A new code should
use the uncondensed system.

### Non-symmetric matrices

For a non-symmetric matrix (e.g. convection-diffusion), give also the block `A_FT` (faces x cells):
`SetupNonsymmetric(A_TT, A_TF, A_FT, A_FF, ...)` then `Solve(b_T, b_F, x_T, x_F)`, or
`SetupNonsymmetricCondensed(A, A_TT, A_TF, A_FT[, A_FF], ...)` then `Solve(b, x)`. The same algorithm runs on the
matrix itself, with GCR instead of FCG and Gauss-Seidel smoothers that follow the direction of the flow (downwind
order). `SetupNonsymmetric()` keeps a reference to `A_FT`, which each `Solve(b_T, b_F, x_T, x_F)` reads
(`b_F - A_FT A_TT^{-1} b_T`): keep it alive and unchanged until the next setup.

### Several processes (MPI)

Compiled with MPI (`-DENABLE_MPI=ON`; `FHHOS4_MPI` is then 1), a solver works on the processes of its communicator,
`Comm`, which must be set (no default):

- `MPI_COMM_SELF`: the system of each process, solved by it alone, as without MPI;
- a communicator of several processes: one system distributed by rows over them. Every setup and solve is then
  collective: called by all the processes of `Comm`, with the same parameters.

The distributed system: each process owns cells and faces, numbered contiguously in rank order (the cells of the
process of rank 0 first, then those of rank 1...; the same for the faces), and gives the rows of its own:

```cpp
solver.Comm = MPI_COMM_WORLD; // to set
// A_TT, A_TF: the rows of the cells of the process; A_FF: the rows of its faces, in full (both triangles).
// The columns of A_TF and A_FF: the faces of the process, then its ghosts (the faces of other processes that its rows
// reference), whose global numbers are faceGhosts, in their local order.
solver.SetupSymmetric(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne, fhhos4::GhostNumbers(nGhosts, faceGhosts));
fhhos4::Result result = solver.Solve(b_T, b_F, x_T, x_F); // the vectors of the cells and faces of the process
```

- **The faces of a process**: its own faces first, in the order of the rows of `A_FF`, then its ghosts, in any order:
  one numbering for the columns of `A_TF` and `A_FF`, whose ghosts are all the faces of other processes in these rows
  (those of its cells, and those coupled to its faces in `A_FF` through the cells of other processes). This is the
  local numbering of DOLFINx's index maps, and of the columns of PETSc's and hypre's distributed matrices.
- **Each face belongs to a process that owns one of its cells**, e.g. the one of the smaller rank; a face has at most
  two cells.
- The vectors (`b_T`, `b_F`, `x_T`, `x_F`, the interpolations of 1): the cells and faces of the process. `Result`: of
  the whole system.
- The algorithm: each process aggregates its own cells (the faces between processes are kept); the coarsest level is
  gathered on every process, which solves it. The number of iterations may change slightly with the number of
  processes.
- Available on several processes: the uncondensed symmetric system (`SetupSymmetric()`), with all the parameters but
  `Krylov::BiCGSTAB`; the Gauss-Seidel smoothers are Gauss-Seidel in each process (the hybrid ones: in each thread),
  Jacobi between the processes. Not yet: `SetupNonsymmetric()`, the condensed interface (`Setup...Condensed()`),
  `CondensedMatrix()`, the Python module.
- The errors of the inputs are raised on all the processes.

### Python

The Python module `fhhos4` has the same interface, with the same names, on NumPy arrays read in place, as the library
reads its inputs: never copied nor converted. It needs Python (>= 3.9), [nanobind](https://nanobind.readthedocs.io) and
NumPy (SciPy is optional), and is built with the library (not installed yet: use it from the build tree):

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_PROGRAM=OFF -DENABLE_PYTHON=ON
cmake --build build
export PYTHONPATH=$PWD/build/python    # the module fhhos4
ctest --test-dir build                 # its tests (they need SciPy)
```

```python
import numpy as np
import fhhos4

solver = fhhos4.Solver()
solver.Dimension = 3   # to set
solver.FaceDegree = k  # to set
solver.CellDegree = k  # to set
solver.SetupSymmetric(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne)  # e.g. scipy.sparse.csr_array, NumPy vectors
x_T, x_F = np.zeros(b_T.size), np.zeros(b_F.size)
result = solver.Solve(b_T, b_F, x_T, x_F)  # x_T, x_F written in place
```

- The parameters are attributes, the enums of `fhhos4.h` classes of the module: `solver.Krylov = fhhos4.Krylov.GCR`.
- The matrices: a `fhhos4.CSRView(rows, cols, rowPointers, columnIndices, values)` of NumPy arrays (or of any
  contiguous array with the buffer protocol), a `fhhos4.CSRMatrix` (the result of the
  helpers), or any object with the attributes `indptr`, `indices`, `data`, `shape` and `format` `"csr"` (SciPy's CSR
  matrices). The indices: int32 or int64, uint64 read as int64; the values: float64. The setups take the index type
  of the library (int32 by default) for both arrays; the other formats raise `TypeError`: convert them with the
  helpers (copies), e.g. `fhhos4.ToCSR(M)` (SciPy CSC matrices included), or view a symmetric CSC matrix with
  `fhhos4.TransposedView(A)`. The CSR format of the library requires sorted columns (SciPy: `sort_indices()`, or
  `sum_duplicates()`, in place).
- The vectors: NumPy vectors of `float64`, read (outputs: written) in place; another dtype raises `TypeError`.
  `Solve(b)` returns `x`.
- The solver keeps the arrays that the solves read alive: those of `A` (setups from the condensed matrix) and of `A_FT`
  (`SetupNonsymmetric()`). `CondensedMatrix()` returns a copy (`fhhos4.CSRMatrix`). The arrays of a `fhhos4.CSRMatrix`
  (`RowPointers`, `ColumnIndices`, `Values`) are read-only NumPy views of the library's memory, which keep it alive.
- The helpers have the same names: `ExtractSymmetricBlocks(M, cells, faces)` with `fhhos4.DoFLists(count, offsets,
  dofs)` (NumPy arrays of the same integer dtype), `Gather()`, `Scatter()`, `ToCSR()`, `TransposedView()`.
- The errors raise `fhhos4.Error`, a `RuntimeError`.
- The setups, solves and helpers release the GIL.

### Inputs

- **The blocks** `A_TT` (cells x cells), `A_TF` (cells x faces), `A_FF` (faces x faces), and `A_FT` (faces x cells) for
  a non-symmetric matrix. The DoFs of each cell (resp. face) are contiguous, cells and faces in any order. In a
  symmetric system, `A_TT` and `A_FF` are read from their lower triangular part. With the condensed interface, `A`
  (both triangles), and `A_FF` is optional.
- **The format of the matrices**: CSR (compressed sparse row), with the index type of the library (`fhhos4::SparseIndex`:
  `int`, or 64 bits with `-DSMALL_INDEX=OFF`) for the row pointers and the columns, the columns strictly increasing in
  each row (sorted, no duplicates; checked). The setups take a `fhhos4::CSRView(rows, cols, rowPointers,
  columnIndices, values)` of your arrays; `fhhos4::SparseMatrix` (Eigen, row-major, with this index type), compressed,
  converts to it. **The matrices are read in place, never copied nor converted**: assemble them in this format. The
  other formats go through the helpers below, which copy them (except `TransposedView()`), so that the copies stay
  visible in your code.
- **The interpolation of the function 1** on the cell and face bases: one coefficient per cell (`cellInterpOfOne`, in
  the order of the blocks of rows of `A_TF`) and per face (`faceInterpOfOne`, blocks of columns). The bases must be
  **hierarchical** (ordered by degree, except with `PCoarsening::ToZero`) **with a constant first function** `phi_0`:
  `1 = c phi_0`, i.e. `c = 1/phi_0`, which may vary from one
  cell or face to another (e.g. L2-orthonormal bases: `c = sqrt(|T|)` on a cell `T`); `c = 1` everywhere if `phi_0 = 1`.
- **The dimension and the degrees**: `Dimension` (2 or 3), `FaceDegree` (k) and `CellDegree` (l) must be set. The setup checks
  the sizes and structure of the inputs against them (full polynomial spaces).

### Other formats: the helpers

| Your matrices | Helper | |
|---|---|---|
| The global matrix of the system (the cells and faces DoFs in your numbering, the Dirichlet DoFs included) | `ExtractSymmetricBlocks(M, cells, faces)`, `ExtractNonsymmetricBlocks(M, cells, faces)`; `Gather()` and `Scatter()` for the vectors | Copies |
| Other index types (64 bits with the 32-bit library, or row pointers and columns of different types) | `ToCSR(M)` | A copy |
| Column-major (Eigen's default, SciPy's CSC) | `ToCSR(M)` (transposed in parallel) | A copy |
| Column-major and symmetric (`A_TT`, `A_FF`, or `A` of a symmetric system) | `TransposedView(A)`: the CSR view of `A^T`, which is `A`, up to the rounding errors of its assembly | No copy |

The copies are `fhhos4::CSRMatrix` objects, owned by the library, which convert to a `CSRView` (valid while the object
lives). From a global matrix (e.g. assembled by a framework in its own numbering):

```cpp
// The DoFs of each cell (resp. face) in the numbering of M: those of the i-th are DoFs[Offsets[i]], ...,
// DoFs[Offsets[i + 1] - 1], in the order of its basis; the DoFs in no list (e.g. of the Dirichlet faces) are dropped
fhhos4::DoFLists<std::int64_t> cells{ nCells, cellOffsets, cellDoFs }, faces{ nFaces, faceOffsets, faceDoFs };
fhhos4::BasicCSRView<std::int64_t, std::int32_t> M(n, n, rowPointers, columns, values); // any index types
fhhos4::HybridBlocks blocks = fhhos4::ExtractSymmetricBlocks(M, cells, faces); // A_TT, A_TF, A_FF: one parallel pass
solver.SetupSymmetric(blocks.A_TT, blocks.A_TF, blocks.A_FF, cellInterpOfOne, faceInterpOfOne);
fhhos4::Gather(b, cells, faces, b_T, b_F); // b: in the numbering of M
solver.Solve(b_T, b_F, x_T, x_F);
fhhos4::Scatter(x_T, x_F, cells, faces, x); // the DoFs of x in no list unchanged (e.g. the Dirichlet values)
```

- `M`: the whole matrix (both triangles), square, in the CSR format (any index types: `std::int32_t` or `std::int64_t`,
  independently).
- The lists: `std::int32_t` or `std::int64_t`; each DoF in one list at most. The order of the lists gives the order of
  the blocks, and that of the DoFs of each cell or face the order of its basis (reorder them, e.g. by degree, here).
- The Dirichlet DoFs, in no list: move their contribution to the right-hand side first (`b - M x_D`).
- `ExtractSymmetricBlocks()` gives the lower triangular parts of `A_TT` and `A_FF` (what `SetupSymmetric()` reads),
  `ExtractNonsymmetricBlocks()` the four blocks in full.

### Parameters

The public members of `fhhos4::Solver`, read by the setup (`Tolerance` and `MaxIterations` by `Solve()`), with the
library's defaults (`library/include/fhhos4.h`, which details each one):

| Parameter | Default | Values |
|---|---|---|
| `Krylov` | `Krylov::Default`: `Krylov::FCG` for a symmetric matrix, `Krylov::GCR` otherwise | `fhhos4::Krylov`: `FCG` (flexible conjugate gradient), `GCR` (GCR(10), any matrix), `BiCGSTAB` (V- or W-cycle only), `None` (the multigrid cycles alone) |
| `Tolerance` | `1e-8` | `Solve()` stops when the relative residual of the condensed system (the true one, `Result::RelativeResidual`) is below it |
| `MaxIterations` | `200` | |
| `Cycle` | `Cycle::K` | `fhhos4::Cycle`: `V`, `W` or `K` (K-cycle: a non-linear preconditioner, which requires FCG or GCR) |
| `PCoarsening` | `PCoarsening::Minus2` | `fhhos4::PCoarsening`, the degrees of the p-levels from `FaceDegree` k down to 0: `Minus2` (k, k-2, ..., 0), `Minus1` (k, k-1, ..., 0), `DivideBy2` (k, k/2, ..., 0), `ToZero` (k, 0: the only p-level keeps the first function of each block, for bases not ordered by degree; usually more iterations at k >= 3) |
| `PreSmoothingIterations`, `PostSmoothingIterations` | `-1`, `-1`: the smoother's, `2` for Chebyshev (its degree), `1` otherwise | `>= 0` (`0`: none) |
| `PreSmoother` | `Smoother::Default`: `BlockChebyshevFourthKind` for a symmetric matrix; for a non-symmetric one, `HybridBlockGaussSeidel` on several threads, `BlockGaussSeidel` on 1 | `fhhos4::Smoother` (one block per face): `BlockGaussSeidel`, `BackwardBlockGaussSeidel`, `SymmetricBlockGaussSeidel` (sequential); `HybridBlockGaussSeidel`, `HybridBackwardBlockGaussSeidel` (their parallel versions); `DampedBlockJacobi` (damping factor 2/3); `BlockChebyshevFirstKind`, `BlockChebyshevFourthKind` (a polynomial in `D^-1 A`, `D` the block diagonal, whose degree is the number of smoothing iterations; parallel, with results that don't depend on the number of threads; symmetric matrices only) |
| `PostSmoother` | `Smoother::Default`: the adjoint of `PreSmoother`, which makes the cycle symmetric: the other direction for Gauss-Seidel (`BackwardBlockGaussSeidel` after `BlockGaussSeidel` and the other way around, the same for the hybrid ones), `PreSmoother` itself if it is symmetric (Chebyshev, Jacobi, `SymmetricBlockGaussSeidel`) | The same as `PreSmoother` |
| `CoarseningFactor` | `3.8` | Ratio of the numbers of cells of two successive h-levels |
| `CoarseMatrixMaxSize` | `1000` | The coarsening stops when the matrix has at most this many rows |
| `CoarseSolver` | `CoarseSolver::Default`: `Cholesky` for a symmetric matrix, `LU` otherwise | `fhhos4::CoarseSolver`: `Cholesky`, `LU` (sparse factorizations) |
| `Threads` | `0`: the OpenMP setting, but at most one thread per physical core unless `OMP_NUM_THREADS` or OpenMP's binding is set (on several processes: the cores divided among the processes of `Comm` that run on the same CPUs) | |
| `Comm` | (MPI only) none: to set | `MPI_COMM_SELF`, or the communicator of the processes of a distributed system |
| `Verbosity` | `0`: silent | `1`: one line per setup and per solve; `2`: all the messages, with the residual at each iteration |

## Constraints

- **Symmetric positive definite** condensed system (`SetupSymmetric...()`), or non-symmetric with a positive definite
  symmetric part (`SetupNonsymmetric...()`): not the pure Neumann problem (singular).
- **Eigen 3.4 or 5.x** in your code, for the header `<fhhos4>`. No Eigen object crosses the boundary of the library:
  its compiled functions take arrays (views), which the header builds from your Eigen matrices and vectors, with your
  version of Eigen and your flags (any version, any alignment settings).
- **The matrices read in place by the solves**: `A` of the setups from the condensed matrix, `A_FT` of
  `SetupNonsymmetric()`. Keep them alive and unchanged until the next setup or the destruction of the solver. The
  others are only read by the setup.
- **Errors**: `fhhos4::Error` (a `std::runtime_error`). `Solve()` does not throw if the tolerance is not reached: see
  `Result::Converged`.
- **Not thread-safe**: one `fhhos4::Solver` object must not be used by several threads at the same time. The library
  parallelizes its own work with OpenMP.
- **MPI**: the calling code initializes MPI (the library calls MPI from the main thread only: `MPI_THREAD_FUNNELED`),
  and `Comm` stays valid while the solver is used.

## License

LGPL-3.0-or-later, like the rest of fhhos4 (see [README.md](README.md)). The library only depends on Eigen (MPL-2.0),
not on the GPL dependencies of the program (GMSH, CGAL): codes under other licenses can link it, under the conditions
of section 4 of the LGPL (in short: keep the license notices, and let the user replace the shared library).

## Version

- `find_package(fhhos4 1.4 REQUIRED)` accepts the versions 1.x with x >= 4: the API stays source-compatible within a
  major version.
- Binary compatibility only within a minor version (`libfhhos4.so.1.4`): adding a parameter to `fhhos4::Solver`
  changes its size. Rebuild your code when fhhos4 changes minor version.
- `FHHOS4_VERSION` (`"1.4.0"`), `FHHOS4_VERSION_MAJOR`, `_MINOR`, `_PATCH`: the version of the header (`<fhhos4>`).
- `fhhos4::Version()`: the version of the library, followed by the commit (`git describe`) when it was built from a git
  clone, e.g. `1.4.0 (v1.4.0-12-gabc1234-dirty)`. Give it in bug reports.
