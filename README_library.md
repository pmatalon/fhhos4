# The fhhos4 multigrid solver as an external library

1. [Installation](#1-installation)
2. [Applicability](#2-applicability)
3. [The algorithm in a nutshell](#3-the-algorithm-in-a-nutshell)
4. [How to use](#3-how-to-use)
   - [Symmetric systems](#symmetric-systems)
   - [Non-symmetric systems](#non-symmetric-systems)
   - [Use as a preconditioner](#use-as-a-preconditioner)
   - [With MPI](#with-mpi)
   - [The Python interface](#the-python-interface)
   - [Inputs](#inputs)
   - [Other matrix formats](#other-matrix-formats)
   - [Parameters](#parameters)
   - [Constraints](#constraints)
5. [License](#license)
6. [Version](#version)

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

To update to a later release `vX.Y.Z`, run `git fetch --tags && git checkout vX.Y.Z` in the clone, then build and install again. `git pull` does nothing in this clone: it follows the tag only.

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

# 2. Applicability

- [x] Linear systems:
   - [x] Symmetric positive definite (e.g. diffusion)
   - [x] Non-symmetric, indefinite systems (e.g. convection-diffusion, both diffusion- and convection-dominated)
   - [x] Singular, symmetric systems
      - [x] with translation kernel (e.g. pure Neumann diffusion)
      - [ ] user-defined kernel vectors (e.g. for elasticity problems)
- [x] Vector-valued unknowns
- [x] `P` and `Q` polynomial spaces in the cells and faces
- [x] Arbitrary polynomial degree:
   - [x] Hierarchical bases **ordered by degree**
   - [x] If not ordered by degree: the degree 0 must at least be first. Then OK up to degree 2.
   - [x] Helper to reorder the bases
   - [ ] Non-hierarchical bases not managed

# 3. The algorithm in a nutshell

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

# 3. How to use

Your code assembles the uncondensed system, i.e.
- the matrix blocks `A_TT, A_TF, A_FF` for a symmetric system, and `A_FT` for a non-symmetric one,
- the two blocks of the right-hand side `(b_T, b_F)`, without the DoFs of the Dirichlet faces (moved to the right-hand side).

That's it, no static condensation.
The library condenses the system internally and solves it for both the
cell and face unknowns `(x_T, x_F)`.

Krylov methods are included: `fhhos4` is meant to be used directly as a solver. It can also be used as a preconditioner for an outer Krylov method.

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

The blocks are Eigen matrices, read
in place, in the format of the library: `fhhos4::SparseMatrix` (row-major, compressed), or an Eigen map of your own
arrays (see [Inputs](#inputs); the other Eigen matrices convert by Eigen, a copy).

### Non-symmetric systems

For a non-symmetric matrix (e.g. convection-diffusion), give also the block `A_FT` by using the `SetupNonSym` function.
```cpp
solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne);
```
The same algorithm runs, with GCR instead of FCG, and Gauss-Seidel smoothers that follow the direction of the flow (downwind
order). The default smoother, `L1HybridBlockGaussSeidel`, converges whatever the numbers of threads and processes. For
convection-dominated problems, `HybridBlockGaussSeidel` usually takes fewer iterations; keep the default if the
iterations stall, e.g. in 3D at degree 2 or more, or on several processes:
```cpp
solver.PreSmoother = fhhos4::Smoother::HybridBlockGaussSeidel;
```

### Use as a preconditioner

The library includes everything to use it as a solver, Krylov method included. 
However, it can also precondition an external Krylov
method working on the condensed system:

```cpp
solver.Krylov = fhhos4::Krylov::None;                 // disable Krylov method in fhhos4
solver.MaxIterations = 1;                             // one cycle per application
solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne);
Eigen::VectorXd b(b_F.size());
solver.CondenseRHS(b_T, b_F, b);                      // b: the condensed right-hand side

// Your Krylov method on A x_F = b
Eigen::VectorXd x_F = Eigen::VectorXd::Zero(b.size()), r = b, z(b.size()), p, Ap(b.size());
while (r.norm() > tolerance * b.norm())
{
    z.setZero();
    solver.Solve(r, z); // the preconditioner, z = M^{-1} r: one cycle
    p = ...;            // the new search direction, from z
    solver.ApplyCondensedMatrix(p, Ap); // the matrix-vector product, Ap = A p
    ...                 // the updates of x_F and r
}

Eigen::VectorXd x_T(b_T.size());
solver.SolveCellUnknowns(b_T, x_F, x_T); // use fhhos4 to solve the cell unknowns at the end
```

On several processes, the same calls with the vectors of the faces of each process (the dot products of your Krylov
method summed over the processes). On one process, `CondensedMatrix()` gives the matrix itself (an Eigen map).

The default K-cycle is non-linear, which means that your Krylov method **must be flexible** (FCG, FGMRES, GCR).
V- and W-cycles work with any Krylov method (CG, GMRES), but are not recommended for performance reasons.

### With MPI

Compiled with MPI (`-DENABLE_MPI=ON`), the solver works on the processes of its communicator,
`Comm`:

- not set (`MPI_COMM_NULL`, the default) or `MPI_COMM_SELF`: the system of the calling process, solved by it alone, as
  without MPI. Not set, the solver makes no MPI call: a code that doesn't use MPI needs no change, MPI not even
  initialized (but for `CoarseSolver::MUMPS`, which needs it).
- a communicator of several processes: the system is distributed by rows.

The distributed system: each process owns cells and faces, numbered contiguously in rank order (the cells of the
process of rank 0 first, then those of rank 1...; the same for the faces), and gives the rows of its own:

```cpp
solver.Comm = MPI_COMM_WORLD;
// A_TT, A_TF: the rows of the cells of the process; A_FF: the rows of its faces, in full (both triangles). The columns
// of A_TF and A_FF: the nOwn faces of the process, then its ghosts (the faces of other processes that its rows
// reference). faceGhostGlobalNumbers[i]: the global number of the ghost of local number nOwn + i (faces, not DoFs).
// E.g. rank 0 owns the faces 0 to 9 and its rows reference the faces 15 and 12 of rank 1, its local faces 10 and 11:
// faceGhostGlobalNumbers = {15, 12}.
std::vector<std::int64_t> faceGhostGlobalNumbers = ...; // or a std::span, an Eigen vector..., GhostNumbers(n, pointer)
solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne, faceGhostGlobalNumbers);
fhhos4::Result result = solver.Solve(b_T, b_F, x_T, x_F); // the vectors of the cells and faces of the process
```

- **`A_FF` in full** on several processes (both triangles; a triangular one is refused).
- **The faces of a process**: its own faces first, in the order of the rows of `A_FF`, then its ghosts, in any order:
  one numbering for the columns of `A_TF` and `A_FF`, whose ghosts are all the faces of other processes in these rows
  (those of its cells, and those coupled to its faces in `A_FF` through the cells of other processes). This is the
  local numbering of DOLFINx's index maps, and of the columns of PETSc's and hypre's distributed matrices.
- **Each face belongs to a process that owns one of its cells**, e.g. the one of the smaller rank; a face has at most
  two cells. A process may own no cell nor face (empty matrices and vectors): it still calls the setups and solves.
- **Non-symmetric**: `SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne, faceGhostGlobalNumbers,
  cellGhostGlobalNumbers)`, with `A_TT` and `A_FF` in full; `A_FT`: the rows of the faces of the process, its columns
  the cells of the process (in the order of the rows of `A_TT`), then its ghost cells, in any order (the cells of other
  processes in these rows), whose global numbers are `cellGhostGlobalNumbers` (numbered as the faces: those of rank 0
  first...). `A_FT` is the only block whose columns reach the cells of other processes: `SetupSym()` has no cell
  ghosts. On each process, the Gauss-Seidel smoothers follow the downwind order of its own faces.
- The vectors (`b_T`, `b_F`, `x_T`, `x_F`, the interpolations of 1): the cells and faces of the process. `Result`: of
  the whole system.
- The algorithm: each process aggregates its own cells (the faces between processes are kept); the coarsest level is
  factorized and solved by MUMPS on all the processes, distributed. The number of iterations may change slightly with
  the number of processes.
- Available on several processes: the uncondensed system (`SetupSym()`, `SetupNonSym()`), with all the
  parameters; the Gauss-Seidel smoothers are Gauss-Seidel in each process (the hybrid ones: in each thread), Jacobi
  between the processes, which can diverge (e.g. in 3D at degree 2 or more), except in the l1 ones (the default).
  `CondensedMatrix()`: one process only (its products: `ApplyCondensedMatrix()`).
- **Global columns**: if your blocks have the global numbers of the faces (and cells) as columns, `LocalColumnsSym()`
  and `LocalColumnsNonSym()` (copies, no communication) number them locally and give the ghosts:
  ```cpp
  // A_TF, A_FF: the rows of the process; the DoF k of the face f is the column f * faceBlockSize + k
  fhhos4::DistributedBlocks blocks = fhhos4::LocalColumnsSym(A_TF, A_FF, firstFace, faceBlockSize);
  solver.SetupSym(A_TT, blocks.A_TF, blocks.A_FF, cellInterpOfOne, faceInterpOfOne, blocks.FaceGhostGlobalNumbers);
  // Non-symmetric: LocalColumnsNonSym(A_TF, A_FT, A_FF, firstFace, firstCell, faceBlockSize, cellBlockSize), whose
  // CellGhostGlobalNumbers are the cellGhostGlobalNumbers of SetupNonSym()
  ```
  `firstFace`: the global number of the first face of the process (`firstCell`, of its first cell); the ghosts by
  increasing global number.
- **Threads**: with several threads per process, give MUMPS's BLAS one thread (e.g. `BLIS_NUM_THREADS=1` with
  conda-forge's MUMPS, which uses BLIS; `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1` for others): a threaded BLAS
  (BLIS starts its threads at each call) makes the many small products of the coarsest solves much slower. MUMPS
  itself runs on one thread (its matrix is small). The same for `CoarseSolver::MUMPS` on one process.
- In Python: see [Python](#python).
- The errors of the inputs are raised on all the processes.

### The Python interface

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
  dofs)` (NumPy arrays of the same integer dtype), `Gather()`, `Scatter()`, `ToCSR()`, `TransposedView()`,
  `LocalColumnsSym(A_TF, A_FF, firstFace, faceBlockSize)` and `LocalColumnsNonSym()` (their `FaceGhostGlobalNumbers` and
  `CellGhostGlobalNumbers`: read-only `int64` NumPy arrays, the `faceGhostGlobalNumbers` and `cellGhostGlobalNumbers` of
  the setups).
- The errors raise `fhhos4.Error`, a `RuntimeError`.
- The setups, solves and helpers release the GIL.
- Compiled with MPI (`-DENABLE_MPI=ON` too; `fhhos4.WithMPI` is then `True`), `solver.Comm` not set: one process, as
  in C++ (see [Several processes](#several-processes-mpi)). To set it, the module needs
  [mpi4py](https://mpi4py.readthedocs.io), built with the same MPI as the library, and imported first (it initializes
  MPI). `solver.Comm`: an mpi4py communicator (`MPI.COMM_SELF`, `MPI.COMM_WORLD`, or its Fortran handle
  `comm.py2f()`). On several processes, `SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne,
  faceGhostGlobalNumbers)`: `faceGhostGlobalNumbers`, the global numbers of the ghost faces, a NumPy array of `int64`
  read in place (`None` by default: no ghosts); `SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne,
  faceGhostGlobalNumbers, cellGhostGlobalNumbers)`, `cellGhostGlobalNumbers` likewise.

  ```python
  from mpi4py import MPI  # first: initializes MPI
  import fhhos4

  solver = fhhos4.Solver()
  solver.Comm = MPI.COMM_WORLD
  ...
  solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne, faceGhostGlobalNumbers)  # the rows of the process
  result = solver.Solve(b_T, b_F, x_T, x_F)  # the vectors of the cells and faces of the process
  ```

### Inputs

- **The blocks** `A_TT` (cells x cells), `A_TF` (cells x faces), `A_FF` (faces x faces), and `A_FT` (faces x cells) for
  a non-symmetric matrix. The DoFs of each cell (resp. face) are contiguous, cells and faces in any order. In a
  symmetric system on one process, `A_TT` and `A_FF` are read from their lower triangular part (lower-only or full);
  otherwise in full (`A_FF` on several processes, all the blocks of a non-symmetric system). A triangular block where
  both triangles are read (or an upper triangular one where the lower is) is refused.
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

### Other matrix formats

`fhhos4` provides the following helper functions to convert your matrices to the format `fhhos4` expects.

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
| `PreSmoother` | `Smoother::Default`: `BlockChebyshevFourthKind` if symmetric, else `L1HybridBlockGaussSeidel` | One block per face: <ul><li>`BlockGaussSeidel`, `BackwardBlockGaussSeidel`, `SymmetricBlockGaussSeidel`: sequential</li><li>`HybridBlockGaussSeidel`, `HybridBackwardBlockGaussSeidel`: parallel</li><li>`L1HybridBlockGaussSeidel`, `L1HybridBackwardBlockGaussSeidel`: parallel, l1 Gauss-Seidel (the diagonal of the rows coupled to other threads or processes increased by `L1Theta` times the sum of the \|a_ij\| of these couplings): convergent whatever the numbers of threads and processes</li><li>`DampedBlockJacobi`: damping 2/3</li><li>`BlockChebyshevFirstKind`, `BlockChebyshevFourthKind`: parallel, symmetric matrices only</li></ul> |
| `PostSmoother` | `Smoother::Default`: the adjoint of `PreSmoother` (symmetric cycle) | Same as `PreSmoother` |
| `L1Theta` | `0.5` | `>= 0`: the factor of the l1 terms of the `L1Hybrid...` smoothers (the sum of the \|a_ij\| of the couplings to other threads or processes added to the diagonal). From `0.5`, convergent for a symmetric matrix whatever the partition; smaller values may save iterations where the hybrid sweeps converge (`0`: `HybridBlockGaussSeidel`) |
| `Prolongation` | `Prolongation::HarmonicExtension` | Values of the faces interior to an aggregate: <ul><li>`HarmonicExtension`: solve their rows of the condensed matrix (shorter setup)</li><li>`Reconstruction`: the article's, from the coarse cell, then one damped Jacobi iteration</li></ul> |
| `CoarseningFactor` | `3.8` | Ratio of the numbers of cells of two successive h-levels |
| `CoarseMatrixMaxSize` | `1000` | Maximum number of rows of the coarsest matrix |
| `CoarseSolver` | `Cholesky` if symmetric, `LU` if non-symmetric | <ul><li>`Cholesky`, `LU`: Eigen's on one process, MUMPS's on several</li><li>`MUMPS`: MUMPS's on one process too</li></ul> |
| `MixedPrecision` | `true` | Levels stored and cycled in single precision (half the memory, faster); the Krylov method, residuals and solution stay in double. Chebyshev and Gauss-Seidel smoothers (all but `DampedBlockJacobi`); otherwise ignored |
| `BlockCSRStorage` | `false` | Level operators stored by blocks of faces when `FaceDegree >= 1` (faster smoothing). Always so in mixed precision; in double, costs a copy of the condensed matrix. Chebyshev and Gauss-Seidel smoothers |
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

# License

LGPL-3.0-or-later, like the rest of fhhos4 (see [README.md](README.md)). The library only depends on Eigen (MPL-2.0),
and, compiled with MPI or `ENABLE_MUMPS`, on MPI and MUMPS (CeCILL-C, compatible with the LGPL; MUMPS's own
dependencies depend on its build: conda-forge's `mumps-mpi` brings ParMETIS, whose license forbids commercial use), not
on the GPL dependencies of the program (GMSH, CGAL): codes under other licenses can link it, under the conditions
of section 4 of the LGPL (in short: keep the license notices, and let the user replace the shared library).

# Version

- `find_package(fhhos4 1.5 REQUIRED)` accepts the versions 1.x with x >= 5.
- Binary compatibility only within a minor version (`libfhhos4.so.1.5`): adding a parameter to `fhhos4::Solver`
  changes its size. Rebuild your code when fhhos4 changes minor version.
- `FHHOS4_VERSION` (`"1.5.0"`), `FHHOS4_VERSION_MAJOR`, `_MINOR`, `_PATCH`: the version of the header (`<fhhos4>`).
- `fhhos4::Version()`: the version of the library, followed by the commit (`git describe`) when it was built from a git
  clone, e.g. `1.5.0 (v1.5.0-12-gabc1234-dirty)`. Give it in bug reports.
