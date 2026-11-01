# The fhhos4 solver as an external library

1. [Installation](#1-installation)
   - [Install from source](#install-from-source)
   - [Conda and Spack packages](#conda-and-spack-packages)
   - [Link](#link)
2. [The algorithm in a nutshell](#2-the-algorithm-in-a-nutshell)
   - [Applicability](#applicability)
3. [How to use](#3-how-to-use)
   - [Symmetric systems](#symmetric-systems)
   - [Non-symmetric systems](#non-symmetric-systems)
   - [As a preconditioner](#as-a-preconditioner)
   - [Several processes (MPI)](#several-processes-mpi)
   - [Python](#python)
   - [Inputs](#inputs)
   - [Other matrix formats: provided helpers](#other-matrix-formats-provided-helpers)
   - [Parameters](#parameters)
   - [Constraints](#constraints)
   - [License](#license)
   - [Version](#version)

# 1. Installation

## Install from source

For use in shared memory only, you'll need
- a C++17 compiler,
- CMake (>= 3.20),
- Ninja or Make,
- [Eigen](https://eigen.tuxfamily.org) 3.4 or 5.x.

For distributed memory, you'll also need
- MPI
- [MUMPS](https://mumps-solver.org)  (e.g. conda-forge's `mumps-mpi`)

```bash
git clone --branch v1.4.0 --depth 1 https://github.com/pmatalon/fhhos4.git && cd fhhos4
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_PROGRAM=OFF
cmake --build build
cmake --install build --prefix <prefix>    # lib/, include/ and the CMake package
```

`v1.4.0` is the latest release. Between releases, this README file on `master` describes the development version; the release's is [README_library.md at v1.4.0](https://github.com/pmatalon/fhhos4/blob/v1.4.0/README_library.md).

Without installing, the library is `build/lib/libfhhos4.so`.

Build options (`cmake -D<OPTION>=<value>`):

| Option | Default | Effect |
|---|---|---|
| `ENABLE_OPENMP` | `ON` | Parallel (OpenMP) |
| `ENABLE_NATIVE_ARCH` | `OFF` | `-march=native`: optimized for the CPU of the build machine, may not run on others |
| `SMALL_INDEX` | `ON` | `int` indices in `fhhos4::SparseMatrix` (`intmax_t` otherwise) |
| `ENABLE_PYTHON` | `OFF` | The Python module, see [Python](#python) |
| `ENABLE_MPI` | `OFF` | The solver on several processes, see [Several processes](#several-processes-mpi). Needs MPI and the MPI version of [MUMPS](https://mumps-solver.org); your code is then compiled and linked with MPI too, and the Python module needs mpi4py |
| `ENABLE_MUMPS` | `OFF` | The coarse solver `CoarseSolver::MUMPS` on one process; always on with `ENABLE_MPI`. Without MPI, needs MUMPS's sequential version, e.g. conda-forge's `mumps-seq` |
| `MUMPS_INCLUDE_DIR`, `MUMPS_LIBRARY` | found | MUMPS is found from its header `dmumps_c.h` and its library `libdmumps`, or given by `-DMUMPS_INCLUDE_DIR=<dir> -DMUMPS_LIBRARY=<path>/libdmumps.so` |

## Conda and Spack packages

Coming soon. Currently waiting for review by the Conda and Spack maintainers.

## Link

```cmake
find_package(fhhos4 REQUIRED)
target_link_libraries(my_code PRIVATE fhhos4::fhhos4)
```

and `-DCMAKE_PREFIX_PATH=<prefix>` (installation), or `-Dfhhos4_DIR=<fhhos4>/build` (build tree). CMake puts the
directory of the library in the RPATH of the executable. A library compiled with MPI brings MPI along
(`MPI::MPI_CXX`).

# 2. The algorithm in a nutshell

The shared library `fhhos4` offers an algebraic multigrid method for the solution of linear systems 
arising from hybrid discretizations (HHO, HDG...) of diffusion and convection-diffusion problems, scalar or
vector-valued, of the form
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
Boundary conditions: Dirichlet (their faces out of the system), Neumann and Robin (their faces in the system), mixed;
the pure Neumann problem (singular: `SingularMatrix`, symmetric systems).
Vector-valued unknowns (`VectorComponents`), e.g. the vector Laplacian: the coarse levels keep the constants of each
component (the translations), not the rotations, which elasticity problems need on their coarse levels (not validated
there).
The polynomial bases used for the discretization **must be hierarchical** and ordered by degree for the p-multigrid part to work.
The polynomial spaces are P^k (simplices, polytopes), or Q^k (quadrilaterals, hexahedra: `PolynomialSpace::Q`), in the
cells and on the faces.
Bases not ordered by degree, but whose first function is constant, work with the p-coarsening straight to the lowest
order (`PCoarsening::ToZero`), or better, reordered by degree: e.g. in the DoF lists of `ExtractBlocksSym()`
(see [Other matrix formats](#other-matrix-formats-provided-helpers)), at no extra cost.

# 3. How to use

Pass the blocks of the uncondensed matrix, `A_TT, A_TF, A_FF` for a symmetric system (and `A_FT` for a non-symmetric
one), and the two blocks of the right-hand side `b_T, b_F`. The library condenses the system and solves it for both the
cell and face unknowns `x_T, x_F`: your code has nothing to do with static condensation, nor with the recovery of the
cell unknowns.

Symmetric positive definite systems, and non-symmetric ones with a positive definite symmetric part
(convection-diffusion). Krylov method included (FCG, GCR). Parallel (OpenMP), and on several processes (MPI,
see [Several processes](#several-processes-mpi)).

### Symmetric systems

```cpp
#include <fhhos4>

fhhos4::Solver solver;
solver.Dimension = 3;  // space dimension
solver.FaceDegree = k; // polynomial degree in the faces
solver.CellDegree = k; // polynomial degree in the cells
// solver.PolynomialSpace = fhhos4::PolynomialSpace::Q; // bases of Q^k instead of P^k
// solver.VectorComponents = 3; // a vector-valued unknown of 3 components (see Inputs)
solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne);
Eigen::VectorXd x_T(b_T.size()), x_F = Eigen::VectorXd::Zero(b_F.size());
fhhos4::Result result = solver.Solve(b_T, b_F, x_T, x_F); // flexible conjugate gradient preconditioned by the AMG
```

Your code assembles the uncondensed system, without the DoFs of the Dirichlet faces (moved to the right-hand side:
`b_T -= A_TD x_D`, `b_F -= A_FD x_D`), and nothing else: no static condensation. The blocks are Eigen matrices, read
in place, in the format of the library: `fhhos4::SparseMatrix` (row-major, compressed), or an Eigen map of your own
arrays (see [Inputs](#inputs); the other Eigen matrices convert by Eigen, a copy). The library condenses the system,
solves the condensed system and computes the cell unknowns, `x_T = A_TT^{-1} (b_T - A_TF x_F)`. The condensed matrix is
available for other uses: `CondensedMatrix()` (an Eigen map), and `Solve(b, x)` solves the condensed system for a
right-hand side `b` condensed by your code or by `CondenseRHS()` (see [As a preconditioner](#as-a-preconditioner)). The setup builds the levels once: call `Solve()` for each right-hand side.
The blocks may be freed after the setup.

### Non-symmetric systems

For a non-symmetric matrix (e.g. convection-diffusion), give also the block `A_FT` by using the `SetupNonSym` function.
```cpp
solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne);
```
The same algorithm runs, with GCR instead of FCG, and Gauss-Seidel smoothers that follow the direction of the flow (downwind
order).

### As a preconditioner

The library includes everything to use it as a solver, Krylov method included. 
However, it can also precondition an external Krylov
method working on the condensed system:

```cpp
solver.Krylov = fhhos4::Krylov::None;                 // disable Krylov method in fhhos4
solver.MaxIterations = 1;                             // one cycle per application
solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne);
fhhos4::SparseMatrixMap A = solver.CondensedMatrix(); // get the condensed matrix computed by fhhos4
Eigen::VectorXd b(b_F.size());
solver.CondenseRHS(b_T, b_F, b);                      // b: the condensed right-hand side

// Your Krylov method on A x_F = b
Eigen::VectorXd x_F = Eigen::VectorXd::Zero(b.size()), r = b, z(b.size()), p, Ap;
while (r.norm() > tolerance * b.norm())
{
    z.setZero();
    solver.Solve(r, z); // the preconditioner, z = M^{-1} r: one cycle
    p = ...;            // the new search direction, from z
    Ap = A * p;         // the matrix-vector product (Eigen)
    ...                 // the updates of x_F and r
}

Eigen::VectorXd x_T(b_T.size());
solver.SolveCellUnknowns(b_T, x_F, x_T); // use fhhos4 to solve the cell unknowns at the end
```

The default K-cycle is non-linear, which means that your Krylov method **must be flexible** (FCG, FGMRES, GCR).
V- and W-cycles work with any Krylov method (CG, GMRES), but are not recommended for performance reasons.

### Several processes (MPI)

Compiled with MPI (`-DENABLE_MPI=ON`; `FHHOS4_MPI` is then 1), a solver works on the processes of its communicator,
`Comm`:

- not set (`MPI_COMM_NULL`, the default) or `MPI_COMM_SELF`: the system of the calling process, solved by it alone, as
  without MPI. Not set, the solver makes no MPI call: a code that doesn't use MPI needs no change, MPI not even
  initialized (but for `CoarseSolver::MUMPS`, which needs it). Forgotten for a distributed system, the setup refuses
  its inputs (ghosts, `A_FF` with more columns than rows) and says to set `Comm`;
- a communicator of several processes: one system distributed by rows over them. Every setup and solve is then
  collective: called by all the processes of `Comm`, with the same parameters.

The distributed system: each process owns cells and faces, numbered contiguously in rank order (the cells of the
process of rank 0 first, then those of rank 1...; the same for the faces), and gives the rows of its own:

```cpp
solver.Comm = MPI_COMM_WORLD;
// A_TT, A_TF: the rows of the cells of the process; A_FF: the rows of its faces, in full (both triangles).
// The columns of A_TF and A_FF: the faces of the process, then its ghosts (the faces of other processes that its rows
// reference), whose global numbers are faceGhosts, in their local order.
solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne, fhhos4::GhostNumbers(nGhosts, faceGhosts));
fhhos4::Result result = solver.Solve(b_T, b_F, x_T, x_F); // the vectors of the cells and faces of the process
```

- **The faces of a process**: its own faces first, in the order of the rows of `A_FF`, then its ghosts, in any order:
  one numbering for the columns of `A_TF` and `A_FF`, whose ghosts are all the faces of other processes in these rows
  (those of its cells, and those coupled to its faces in `A_FF` through the cells of other processes). This is the
  local numbering of DOLFINx's index maps, and of the columns of PETSc's and hypre's distributed matrices.
- **Each face belongs to a process that owns one of its cells**, e.g. the one of the smaller rank; a face has at most
  two cells. A process may own no cell nor face (empty matrices and vectors): it still calls the setups and solves.
- **Non-symmetric**: `SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne, faceGhosts,
  cellGhosts)`, with `A_TT` and `A_FF` in full; `A_FT`: the rows of the faces of the process, its columns the cells of
  the process (in the order of the rows of `A_TT`), then its ghost cells, in any order (the cells of other processes
  in these rows), whose global numbers are `cellGhosts` (numbered as the faces: those of rank 0 first...). On each
  process, the Gauss-Seidel smoothers follow the downwind order of its own faces.
- The vectors (`b_T`, `b_F`, `x_T`, `x_F`, the interpolations of 1): the cells and faces of the process. `Result`: of
  the whole system.
- The algorithm: each process aggregates its own cells (the faces between processes are kept); the coarsest level is
  factorized and solved by MUMPS on all the processes, distributed. The number of iterations may change slightly with
  the number of processes.
- Available on several processes: the uncondensed system (`SetupSym()`, `SetupNonSym()`), with all the
  parameters; the Gauss-Seidel smoothers are Gauss-Seidel in each process (the hybrid ones: in each thread), Jacobi
  between the processes. Not yet: `CondensedMatrix()`.
- In Python: see [Python](#python).
- The errors of the inputs are raised on all the processes.

### Python

The Python module `fhhos4` has the same interface, with the same names, on NumPy arrays read in place, as the library
reads its inputs: never copied nor converted (with its own matrix types, Python having no Eigen: `fhhos4.CSRView`,
`fhhos4.CSRMatrix`, below). It needs Python (>= 3.9), [nanobind](https://nanobind.readthedocs.io) and
NumPy (SciPy is optional), and is built with the library (not installed yet: use it from the build tree):

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_PYTHON=ON
cmake --build build
export PYTHONPATH=$PWD/build/python    # the module fhhos4
ctest --test-dir build                 # its tests (they need SciPy)
```

```python
import numpy as np
import fhhos4

solver = fhhos4.Solver()
solver.Dimension = 3   # space dimension
solver.FaceDegree = k  # polynomial degree in the faces
solver.CellDegree = k  # polynomial degree in the cells
# solver.PolynomialSpace = fhhos4.PolynomialSpace.Q  # bases of Q^k instead of P^k
solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne)  # e.g. scipy.sparse.csr_array, NumPy vectors
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
- The solver keeps the arrays that the solves read alive: those of `A_FT` (`SetupNonSym()`). `CondensedMatrix()`
  returns a copy (`fhhos4.CSRMatrix`). The arrays of a `fhhos4.CSRMatrix`
  (`RowPointers`, `ColumnIndices`, `Values`) are read-only NumPy views of the library's memory, which keep it alive.
- The helpers have the same names: `ExtractBlocksSym(M, cells, faces)` with `fhhos4.DoFLists(count, offsets,
  dofs)` (NumPy arrays of the same integer dtype), `Gather()`, `Scatter()`, `ToCSR()`, `TransposedView()`.
- The errors raise `fhhos4.Error`, a `RuntimeError`.
- The setups, solves and helpers release the GIL.
- Compiled with MPI (`-DENABLE_MPI=ON` too; `fhhos4.WithMPI` is then `True`), `solver.Comm` not set: one process, as
  in C++ (see [Several processes](#several-processes-mpi)). To set it, the module needs
  [mpi4py](https://mpi4py.readthedocs.io), built with the same MPI as the library, and imported first (it initializes
  MPI). `solver.Comm`: an mpi4py communicator (`MPI.COMM_SELF`, `MPI.COMM_WORLD`, or its Fortran handle
  `comm.py2f()`). On several processes, `SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne,
  faceGhosts)`: `faceGhosts`, the global numbers of the ghost faces, a NumPy array of `int64` read in place (`None` by
  default: no ghosts); `SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne, faceGhosts,
  cellGhosts)`, `cellGhosts` likewise.

  ```python
  from mpi4py import MPI  # first: initializes MPI
  import fhhos4

  solver = fhhos4.Solver()
  solver.Comm = MPI.COMM_WORLD
  ...
  solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne, faceGhosts)  # the rows of the process
  result = solver.Solve(b_T, b_F, x_T, x_F)  # the vectors of the cells and faces of the process
  ```

### Inputs

- **The blocks** `A_TT` (cells x cells), `A_TF` (cells x faces), `A_FF` (faces x faces), and `A_FT` (faces x cells) for
  a non-symmetric matrix. The DoFs of each cell (resp. face) are contiguous, cells and faces in any order. In a
  symmetric system, `A_TT` and `A_FF` are read from their lower triangular part.
- **The format of the matrices**: Eigen sparse matrices, compressed, row-major, with the index type of the library
  (`fhhos4::SparseIndex`: `int`, or 64 bits with `-DSMALL_INDEX=OFF`): `fhhos4::SparseMatrix`, or an Eigen map of your
  own CSR arrays, `fhhos4::SparseMatrixMap(rows, cols, nonZeros, rowPointers, columnIndices, values)`, the columns
  strictly increasing in each row (sorted, no duplicates; checked). **The matrices are read in place, never copied nor
  converted**: assemble them in this format. Any other Eigen sparse matrix converts by Eigen, in one line, a copy that
  stays visible in your code: `fhhos4::SparseMatrix B = A;` (see below). A wrong format does not compile, with a
  message that says how to convert.
- **The interpolation of the function 1** on the cell and face bases: one coefficient per cell (`cellInterpOfOne`, in
  the order of the blocks of rows of `A_TF`) and per face (`faceInterpOfOne`, blocks of columns). The bases must be
  **hierarchical** (ordered by degree: for each d, the first functions of a block span P^d, or Q^d whose degree is the
  largest partial degree; except with `PCoarsening::ToZero`) **with a constant first function** `phi_0`:
  `1 = c phi_0`, i.e. `c = 1/phi_0`, which may vary from one
  cell or face to another (e.g. L2-orthonormal bases: `c = sqrt(|T|)` on a cell `T`); `c = 1` everywhere if `phi_0 = 1`.
- **The dimension and the degrees**: `Dimension` (2 or 3), `FaceDegree` (k) and `CellDegree` (l) must be set, and
  `PolynomialSpace` for bases of Q^k. The setup checks the sizes and structure of the inputs against them (full
  polynomial spaces: P^l in the cells and P^k on the faces, or Q^l and Q^k; times `VectorComponents`).
- **A vector-valued unknown** of `n` components: `VectorComponents = n`. Each cell and face block holds the scalar basis
  times each unit vector, **interleaved by component** (array of structures): `phi_i e_c` is the DoF `n i + c` of the
  block. The interpolations of 1 stay those of the scalar bases, one coefficient per cell and per face.

### Other matrix formats: provided helpers

| Your matrices | How to convert | Makes a copy? |
|---|---|---|
| Column-major (Eigen's default), or another index type | `fhhos4::SparseMatrix B = A;` (Eigen's conversion) | Yes |
| Column-major and symmetric (`A_TT`, `A_FF`, or `A` of a symmetric system) | `TransposedView(A)`: `A^T` as an Eigen map of the arrays of `A`, which is `A`, up to the rounding errors of its assembly | No |
| Your own CSR arrays, with the index type of the library | `fhhos4::SparseMatrixMap(rows, cols, nonZeros, rowPointers, columnIndices, values)` | No |
| The global matrix of the system (the cells and faces DoFs in your numbering, the Dirichlet DoFs included) | `ExtractBlocksSym(M, cells, faces)`, `ExtractBlocksNonSym(M, cells, faces)`; `Gather()` and `Scatter()` for the vectors | Yes |

The copies are Eigen matrices of your code (`fhhos4::SparseMatrix`). From a global matrix (e.g. assembled by a framework
in its own numbering):

```cpp
// The DoFs of each cell (resp. face) in the numbering of M: those of the i-th are DoFs[Offsets[i]], ...,
// DoFs[Offsets[i + 1] - 1], in the order of its basis; the DoFs in no list (e.g. of the Dirichlet faces) are dropped
fhhos4::DoFLists<std::int64_t> cells{ nCells, cellOffsets, cellDoFs }, faces{ nFaces, faceOffsets, faceDoFs };
Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor, std::int64_t>> M(n, n, nnz, rowPointers, columns, values);
fhhos4::HybridBlocks blocks = fhhos4::ExtractBlocksSym(M, cells, faces); // A_TT, A_TF, A_FF: one parallel pass
solver.SetupSym(blocks.A_TT, blocks.A_TF, blocks.A_FF, cellInterpOfOne, faceInterpOfOne);
fhhos4::Gather(b, cells, faces, b_T, b_F); // b: in the numbering of M
solver.Solve(b_T, b_F, x_T, x_F);
fhhos4::Scatter(x_T, x_F, cells, faces, x); // the DoFs of x in no list unchanged (e.g. the Dirichlet values)
```

- `M`: the whole matrix (both triangles), square, an Eigen row-major matrix (or a map of your arrays) with the index type
  `std::int32_t` or `std::int64_t`.
- The lists: `std::int32_t` or `std::int64_t`; each DoF in one list at most. The order of the lists gives the order of
  the blocks, and that of the DoFs of each cell or face the order of its basis (reorder them, e.g. by degree, here).
- The Dirichlet DoFs, in no list: move their contribution to the right-hand side first (`b - M x_D`).
- `ExtractBlocksSym()` gives the lower triangular parts of `A_TT` and `A_FF` (what `SetupSym()` reads),
  `ExtractBlocksNonSym()` the four blocks in full.

### Parameters

The public members of `fhhos4::Solver`, read by the setup (`Tolerance` and `MaxIterations` by `Solve()`), with the
library's defaults (`library/include/fhhos4.h`, which details each one):

| Parameter | Default | Values |
|---|---|---|
| `PolynomialSpace` | `PolynomialSpace::P` | Space of the cell and face bases: <ul><li>`P`: degree at most k</li><li>`Q`: degree at most k in each variable (quadrilaterals, hexahedra)</li></ul> |
| `VectorComponents` | `1` | `>= 1`: components of the unknown, interleaved in each block (see [Inputs](#inputs)) |
| `SingularMatrix` | `false` | `true`: the kernel is the interpolation of 1 of each component (pure Neumann; translations, not rotations). `SetupSym()` only. The solution is defined up to that kernel |
| `Krylov` | `FCG` if symmetric, `GCR` if non-symmetric | <ul><li>`FCG`: flexible conjugate gradient</li><li>`GCR`: restarted GCR, any matrix</li><li>`None`: multigrid cycles alone (e.g. [as a preconditioner](#as-a-preconditioner))</li></ul> |
| `Restart` | `1` for FCG, `10` for GCR | `>= 1`: search directions stored (two face vectors each). Fewer save memory; GCR may then need more iterations |
| `Tolerance` | `1e-8` | Relative residual of the condensed system at which `Solve()` stops |
| `MaxIterations` | `200` | |
| `Cycle` | `Cycle::K` | <ul><li>`V`</li><li>`W`</li><li>`K`: non-linear: a Krylov method around it must be flexible (FCG, GCR, FGMRES)</li></ul> |
| `PCoarsening` | `PCoarsening::Minus2` | Degrees of the p-levels, from `FaceDegree` k: <ul><li>`Minus2`: k, k-2, ..., 0</li><li>`Minus1`: k, k-1, ..., 0</li><li>`ToZero`: k, 0, for bases not ordered by degree (more iterations from k = 3)</li></ul> |
| `PreSmoothingIterations`, `PostSmoothingIterations` | `2` for Chebyshev (its degree), else `1` | Any value `>= 0` (`0`: none) |
| `PreSmoother` | `Smoother::Default`: `BlockChebyshevFourthKind` if symmetric, else `HybridBlockGaussSeidel` (`BlockGaussSeidel` on 1 thread) | One block per face: <ul><li>`BlockGaussSeidel`, `BackwardBlockGaussSeidel`, `SymmetricBlockGaussSeidel`: sequential</li><li>`HybridBlockGaussSeidel`, `HybridBackwardBlockGaussSeidel`: parallel</li><li>`DampedBlockJacobi`: damping 2/3</li><li>`BlockChebyshevFirstKind`, `BlockChebyshevFourthKind`: parallel, symmetric matrices only</li></ul> |
| `PostSmoother` | `Smoother::Default`: the adjoint of `PreSmoother` (symmetric cycle) | Same as `PreSmoother` |
| `Prolongation` | `Prolongation::HarmonicExtension` | Values of the faces interior to an aggregate: <ul><li>`HarmonicExtension`: solve their rows of the condensed matrix (shorter setup)</li><li>`Reconstruction`: the article's, from the coarse cell, then one damped Jacobi iteration</li></ul> |
| `CoarseningFactor` | `3.8` | Ratio of the numbers of cells of two successive h-levels |
| `CoarseMatrixMaxSize` | `1000` | Maximum number of rows of the coarsest matrix |
| `CoarseSolver` | `Cholesky` if symmetric, `LU` if non-symmetric | <ul><li>`Cholesky`, `LU`: Eigen's on one process, MUMPS's on several</li><li>`MUMPS`: MUMPS's on one process too</li></ul> |
| `MixedPrecision` | `true` | Levels stored and cycled in single precision (half the memory, faster); the Krylov method, residuals and solution stay in double. Symmetric matrix, Chebyshev smoothers only (in CSR on several processes); otherwise ignored |
| `BlockCSRStorage` | `false` | Level operators stored by blocks of faces when `FaceDegree >= 1` (faster smoothing). Always so in mixed precision; in double, costs a copy of the condensed matrix. Chebyshev smoothers, one process only |
| `Threads` | `0`: OpenMP's setting, at most one per physical core unless `OMP_NUM_THREADS` or OpenMP's binding is set (shared among the processes on the same CPUs) | Number of threads |
| `Comm` | (MPI only) `MPI_COMM_NULL`: one process, no MPI call (MPI may be uninitialized, except with `CoarseSolver::MUMPS`) | <ul><li>`MPI_COMM_SELF`: one process</li><li>the communicator of a distributed system</li></ul> |
| `Verbosity` | `0`: silent | <ul><li>`1`: one line per setup and per solve</li><li>`2`: all messages, with the residual at each iteration</li></ul> |

## Constraints

- **Symmetric positive definite** condensed system (`SetupSym()`), or semi-definite with the interpolation of 1
  (of each component) as kernel (`SingularMatrix`: the pure Neumann problem), or non-symmetric with a positive definite symmetric part
  (`SetupNonSym()`). With `SetupSym()`, the blocks
  of `A_TT` (one per cell) must be positive definite: the library factorizes them by Cholesky (LU with
  `SetupNonSym()`).
- **Eigen 3.4 or 5.x** in your code, for the header `<fhhos4>`. No Eigen object crosses the boundary of the library:
  its compiled functions take arrays (views), which the header builds from your Eigen matrices and vectors, with your
  version of Eigen and your flags (any version, any alignment settings).
- **The matrix read in place by the solves**: `A_FT` of `SetupNonSym()`. Keep it alive and unchanged until the
  next setup or the destruction of the solver. The others are only read by the setup.
- **Errors**: `fhhos4::Error` (a `std::runtime_error`). `Solve()` does not throw if the tolerance is not reached: see
  `Result::Converged`.
- **Not thread-safe**: one `fhhos4::Solver` object must not be used by several threads at the same time. The library
  parallelizes its own work with OpenMP.
- **MPI**: the calling code initializes MPI (the library calls MPI from the main thread only: `MPI_THREAD_FUNNELED`),
  and `Comm` stays valid while the solver is used.

## License

LGPL-3.0-or-later, like the rest of fhhos4 (see [README.md](README.md)). The library only depends on Eigen (MPL-2.0),
and, compiled with MPI or `ENABLE_MUMPS`, on MPI and MUMPS (CeCILL-C, compatible with the LGPL; MUMPS's own
dependencies depend on its build: conda-forge's `mumps-mpi` brings ParMETIS, whose license forbids commercial use), not
on the GPL dependencies of the program (GMSH, CGAL): codes under other licenses can link it, under the conditions
of section 4 of the LGPL (in short: keep the license notices, and let the user replace the shared library).

## Version

- `find_package(fhhos4 1.5 REQUIRED)` accepts the versions 1.x with x >= 5.
- Binary compatibility only within a minor version (`libfhhos4.so.1.5`): adding a parameter to `fhhos4::Solver`
  changes its size. Rebuild your code when fhhos4 changes minor version.
- `FHHOS4_VERSION` (`"1.5.0"`), `FHHOS4_VERSION_MAJOR`, `_MINOR`, `_PATCH`: the version of the header (`<fhhos4>`).
- `fhhos4::Version()`: the version of the library, followed by the commit (`git describe`) when it was built from a git
  clone, e.g. `1.5.0 (v1.5.0-12-gabc1234-dirty)`. Give it in bug reports.
