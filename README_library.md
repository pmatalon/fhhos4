# The fhhos4 multigrid solver as an external library

1. [Installation](#1-installation)
2. [Applicability](#2-applicability)
3. [The algorithm in a nutshell](#3-the-algorithm-in-a-nutshell)
4. [How to use](#3-how-to-use)
   - [Symmetric systems](#symmetric-systems)
   - [Inputs](#inputs)
   - [With MPI](#with-mpi)
   - [Non-symmetric systems](#non-symmetric-systems)
   - [Use as a preconditioner](#use-as-a-preconditioner)
   - [The Python interface](#the-python-interface)
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
git clone --branch v1.5.0 --depth 1 https://github.com/pmatalon/fhhos4.git && cd fhhos4
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_PROGRAM=OFF
cmake --build build
cmake --install build --prefix <prefix>    # lib/, include/ and the CMake package
```

`v1.5.0` is the latest release. Between releases, this README file on `master` describes the development version; the release's is [README_library.md at v1.5.0](https://github.com/pmatalon/fhhos4/blob/v1.5.0/README_library.md).

To update to a later release `vX.Y.Z`, run `git fetch --tags && git checkout vX.Y.Z` in the clone, then build and install again. `git pull` does nothing in this clone: it follows the tag only.

Without installing, the library is `build/lib/libfhhos4.so`.

Build options (`cmake -D<OPTION>=<value>`):

| Option | Default | Effect |
|---|---|---|
| `ENABLE_OPENMP` | `ON` | Parallel (OpenMP) |
| `ENABLE_MPI` | `OFF` | The solver on several processes, see [Several processes](#with-mpi). Needs MPI and the MPI version of [MUMPS](https://mumps-solver.org); your code is then compiled and linked with MPI too, and the Python module needs mpi4py |
| `ENABLE_PYTHON` | `OFF` | The Python module, see [Python](#the-python-interface) |
| `ENABLE_MUMPS` | `OFF` | The coarse solver `CoarseSolver::MUMPS` on one process; always on with `ENABLE_MPI`. Without MPI, needs MUMPS's sequential version, e.g. conda-forge's `mumps-seq` |
| `MUMPS_INCLUDE_DIR`, `MUMPS_LIBRARY` | found | MUMPS is found from its header `dmumps_c.h` and its library `libdmumps`, or given by `-DMUMPS_INCLUDE_DIR=<dir> -DMUMPS_LIBRARY=<path>/libdmumps.so` |
| `ENABLE_NATIVE_ARCH` | `OFF` | `-march=native`: optimized for the CPU of the build machine, may not run on others |
| `ENABLE_64BIT_INDEX` | `OFF` | 64-bit indices (`int64_t`) in `fhhos4::SparseMatrix` (`int` otherwise) |

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
   - [ ] p-adapdability not managed (fixed polynomial degrees across the mesh!)

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

### Symmetric systems

```cpp
#include <fhhos4>

fhhos4::Solver solver;
solver.Dimension = 3;  // space dimension
solver.FaceDegree = k; // polynomial degree in the faces
solver.CellDegree = k; // polynomial degree in the cells
// solver.PolynomialSpace = fhhos4::PolynomialSpace::Q; // bases of Q^k instead of P^k
// solver.VectorComponents = 3; // a vector-valued unknown of 3 components (see Inputs)

#if FHHOS4_MPI // 1 if the library is compiled with MPI (fhhos4_version.h)
solver.Comm = MPI_COMM_WORLD;
// The blocks: the rows of the cells and faces of this process, the columns of the faces numbered locally (see With MPI)
std::vector<std::int64_t> faceGhostGlobalNumbers = ...; // or a std::span, an Eigen vector, GhostGlobalNumbers(n, p)
solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne, faceGhostGlobalNumbers);
#else
solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne);
#endif

Eigen::VectorXd x_T(b_T.size()), x_F = Eigen::VectorXd::Zero(b_F.size()); // with MPI: of the process
fhhos4::Result result = solver.Solve(b_T, b_F, x_T, x_F); // flexible conjugate gradient preconditioned by the AMG
```

The blocks are of type `fhhos4::SparseMatrix`, which is an Eigen sparse matrices (row-major, compressed).
It can also an Eigen map of your own arrays (see [Inputs](#inputs)).

### Inputs

- **The blocks** `A_TT`, `A_TF`, `A_FF`, and `A_FT` for
  a non-symmetric matrix. The DoFs of each cell (resp. face) are contiguous, cells and faces in any order. In a
  symmetric system on one process, `A_TT` and `A_FF` are read from their lower triangular part (lower-only or full);
  otherwise in full (`A_FF` on several processes, all the blocks of a non-symmetric system). A triangular block where
  both lower and upper triangular parts are read (or an upper triangular one where the lower is) is refused.
- **The format of the matrices**: Eigen sparse matrices, compressed, row-major, with the index type of the library
  (`fhhos4::SparseIndex`: `int`, or 64 bits with `-DENABLE_64BIT_INDEX=ON`): `fhhos4::SparseMatrix`, or an Eigen map of your
  own CSR arrays, `fhhos4::SparseMatrixMap(rows, cols, nonZeros, rowPointers, columnIndices, values)`, the columns
  strictly increasing in each row (sorted, no duplicates; checked). **The matrices are read in place, never copied nor
  converted**: assemble them in this format. Any other Eigen sparse matrix converts by Eigen, in one line, a copy that
  stays visible in your code: `fhhos4::SparseMatrix B = A;` (see below).
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

### With MPI

If the library is compiled with MPI (`-DENABLE_MPI=ON`), `solver.Comm` sets the processes that solve the system:

- **Not set** (`MPI_COMM_NULL`, the default) or `MPI_COMM_SELF`: each process solves its own system alone, as without
  MPI. If `Comm` is not set, the solver makes no MPI call: a code that doesn't use MPI needs no change, and doesn't even
  have to initialize MPI (unless it uses `CoarseSolver::MUMPS`).
- **A communicator of several processes**: the processes solve one system together, distributed by rows. Every setup
  and solve is collective: all the processes of the communicator call it, with the same parameters.

**Who owns what.** Each cell and each face belongs to one process. A face must belong to a process that owns one of
its cells (a face has one or two cells), e.g. the one of smaller rank. A process may own no cell and no face: it then
gives empty matrices and vectors, and still calls the setups and solves.

#### Expected format

Each process gives `SetupSym()` the following (the `FHHOS4_MPI` branch of the example above):

- `A_TT` and `A_TF`: the rows of its cells. `A_FF`: the rows of its faces, assembled (with the contributions of the
  cells of other processes) and with both lower and upper triangular parts (on one process, the lower triangular part
  alone is enough).
- The columns of `A_TF` and `A_FF` follow a local numbering of the faces: first the faces of the process, in the order
  of the rows of `A_FF`, then its **ghost faces**, in any order. The ghost faces are the faces of other processes that
  appear in these rows: the faces of its cells that other processes own, and the faces coupled to its own faces in
  `A_FF` through the cells of other processes. This is the local numbering of DOLFINx's index maps.
- `faceGhostGlobalNumbers`: the global numbers of its ghost faces (faces, not DoFs), in the order of their local
  numbers. The faces are numbered as PETSc and hypre require for their distributed matrices: those of each process
  contiguously, in rank order.
- The vectors `b_T`, `b_F`, `x_T`, `x_F` and the interpolations of 1: the values on its own cells and faces.

For example, with two processes, rank 0 owns the faces 0 to 9 and rank 1 the faces 10 to 19. If the rows of rank 0
contain the faces 15 and 12 of rank 1, rank 0 numbers them locally 10 and 11 (or 11 and 10) and gives
`faceGhostGlobalNumbers = {15, 12}` (or `{12, 15}`).

PETSc's and hypre's distributed matrices store each block in two parts: `MergeOffDiag()` converts them (see
[Other matrix formats](#other-matrix-formats)).

#### What you get back

`Solve()` writes `x_T` and `x_F` on the cells and faces of the process, and returns the `Result`
of the whole system. An error in the inputs of one process is raised on all the processes, so that none is left
waiting.

#### What changes in the solver

- Each process aggregates its own cells; the faces between processes are kept. The coarsest level is factorized and
  solved by MUMPS, distributed over all the processes. The number of iterations may change slightly with the number
  of processes.
- All the parameters work on several processes. The Gauss-Seidel smoothers are Gauss-Seidel within each process (the
  hybrid ones: within each thread), and Jacobi between the processes, which can diverge (e.g. in 3D at degree 2 or
  more), except the l1 smoothers (the default).
- `CondensedMatrix()` works on one process only; `ApplyCondensedMatrix()` on any number.
- With several threads per process, give MUMPS's BLAS one thread (`BLIS_NUM_THREADS=1` with conda-forge's MUMPS,
  which uses BLIS; `OPENBLAS_NUM_THREADS=1` or `MKL_NUM_THREADS=1` with others): a threaded BLAS makes the many small
  products of the coarsest solves much slower (BLIS starts its threads at each call). MUMPS itself runs on one thread
  (its matrix is small). The same holds for `CoarseSolver::MUMPS` on one process.

In Python: see [Python](#the-python-interface).

### Non-symmetric systems

What changes for a non-symmetric matrix (e.g. convection-diffusion):

- The block `A_FT` too, given to `SetupNonSym()`:
  ```cpp
  solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne);
  ```
- The same algorithm runs, with GCR instead of FCG, and Gauss-Seidel smoothers that follow the direction of the flow
  (downwind order). The default smoother, `L1HybridBlockGaussSeidel`, converges whatever the numbers of threads and
  processes. For convection-dominated problems, `HybridBlockGaussSeidel` usually takes fewer iterations; keep the
  default if the iterations stall, e.g. in 3D at degree 2 or more, or on several processes:
  ```cpp
  solver.PreSmoother = fhhos4::Smoother::HybridBlockGaussSeidel;
  ```
- [With MPI](#with-mpi), each process gives the rows of its faces in `A_FT`. Its columns follow a local numbering of
  the cells, as those of `A_TF` and `A_FF` for the faces: first the cells of the process, in the order of the rows of
  `A_TT`, then its **ghost cells**, in any order (the cells of other processes that appear in these rows). `A_FT` is
  the only block with ghost cells. Their global numbers, `cellGhostGlobalNumbers`, number the cells as
  `faceGhostGlobalNumbers` the faces:
  ```cpp
  solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne, faceGhostGlobalNumbers,
                     cellGhostGlobalNumbers);
  ```
  On each process, the Gauss-Seidel smoothers follow the downwind order of its own faces.

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

### The Python interface

The Python module `fhhos4` has the classes, functions, parameters and enums of the C++ library, with the same names
(`solver.Krylov = fhhos4.Krylov.GCR`), and also reads its inputs in place, as NumPy arrays. It needs Python (>= 3.9),
[nanobind](https://nanobind.readthedocs.io) and NumPy (SciPy is optional), and is built with the library (not installed
yet: use it from the build tree):

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
solver.Dimension = 3
solver.FaceDegree = k
solver.CellDegree = k
solver.SetupSym(A_TT, A_TF, A_FF, cellInterpOfOne, faceInterpOfOne)  # e.g. scipy.sparse.csr_array, NumPy vectors
x_T, x_F = np.zeros(b_T.size), np.zeros(b_F.size)
result = solver.Solve(b_T, b_F, x_T, x_F)  # x_T, x_F written in place
```

What differs from C++:

- **The matrices**: SciPy's CSR matrices (or any object with the attributes `indptr`, `indices`, `data`, `shape` and
  `format` `"csr"`), a `fhhos4.CSRView(rows, cols, rowPointers, columnIndices, values)` of NumPy arrays (or of any
  contiguous array with the buffer protocol), or a `fhhos4.CSRMatrix` (what the helpers return). The values: float64;
  the columns sorted (SciPy: `sort_indices()`, in place). The setups take the index type of the library (int32 by
  default) for both index arrays, the helpers int32 or int64 (uint64 read as int64). Other formats raise `TypeError`:
  convert them with `fhhos4.ToCSR(M)` (a copy, SciPy's CSC matrices included), or view a symmetric CSC matrix with
  `fhhos4.TransposedView(A)`.
- **The vectors**: NumPy vectors of float64, read (outputs: written) in place; another dtype raises `TypeError`. The
  lists of `fhhos4.DoFLists(count, offsets, dofs)`: NumPy arrays of the same integer dtype.
- **Lifetimes**: the solver keeps alive the arrays that its solves read (those of `A_FT`, for `SetupNonSym()`).
  `CondensedMatrix()` returns a copy. The arrays of a `fhhos4.CSRMatrix` (`RowPointers`, `ColumnIndices`, `Values`)
  are read-only NumPy views of the library's memory, which they keep alive.
- **The errors** raise `fhhos4.Error`, a `RuntimeError`. The setups, solves and helpers release the GIL.
- **With MPI** (`-DENABLE_MPI=ON` too; `fhhos4.WithMPI` is then `True`): `solver.Comm` is an
  [mpi4py](https://mpi4py.readthedocs.io) communicator (`MPI.COMM_WORLD`, or its Fortran handle `comm.py2f()`). Import
  mpi4py, built with the same MPI as the library, before `fhhos4` (it initializes MPI). The ghost numbers are NumPy
  arrays of int64, read in place (`None` by default: no ghosts). Without `Comm`, mpi4py is not needed. The
  `offdColumns` of `MergeOffDiag()`: NumPy arrays of int32 or int64; its ghosts: arrays of int64.

### Other matrix formats

`fhhos4` provides the following helper functions to convert your matrices to the format `fhhos4` expects (see
[Inputs](#inputs)).

| Your matrices | How to convert | Makes a copy? |
|---|---|---|
| Column-major (Eigen's default), or another index type | `fhhos4::SparseMatrix B = A;` (Eigen's conversion) | Yes |
| Column-major and symmetric (`A_TT`, `A_FF`, or `A` of a symmetric system) | `TransposedView(A)`: `A^T` as an Eigen map of the arrays of `A`, which is `A`, up to the rounding errors of its assembly | No |
| Your own CSR arrays, with the index type of the library | `fhhos4::SparseMatrixMap(rows, cols, nonZeros, rowPointers, columnIndices, values)` | No |
| The monolithic matrix of the system (the cells and faces DoFs in your numbering, the Dirichlet DoFs included) | `ExtractBlocksSym(M, cells, faces)`, `ExtractBlocksNonSym(M, cells, faces)`; `SplitVector()` and `JoinVector()` for the vectors | Yes |
| With MPI, blocks stored in two parts, as PETSc (MPIAIJ) and hypre (ParCSR) store them: the columns of the faces of the process, and the others with their global numbers | `MergeOffDiag(...)` (below) | Yes (not `A_TT`) |

The copies are Eigen matrices of your code (`fhhos4::SparseMatrix`). From a monolithic matrix (e.g. assembled by a
framework in its own numbering):

```cpp
// The DoFs of each cell (resp. face) in the numbering of M: those of the i-th are DoFs[Offsets[i]], ...,
// DoFs[Offsets[i + 1] - 1], in the order of its basis; the DoFs in no list (e.g. of the Dirichlet faces) are dropped
fhhos4::DoFLists<std::int64_t> cells{ nCells, cellOffsets, cellDoFs }, faces{ nFaces, faceOffsets, faceDoFs };
Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor, std::int64_t>> M(n, n, nnz, rowPointers, columns, values);
fhhos4::HybridBlocks blocks = fhhos4::ExtractBlocksSym(M, cells, faces); // A_TT, A_TF, A_FF: one parallel pass
solver.SetupSym(blocks.A_TT, blocks.A_TF, blocks.A_FF, cellInterpOfOne, faceInterpOfOne);
fhhos4::SplitVector(b, cells, faces, b_T, b_F); // b: in the numbering of M
solver.Solve(b_T, b_F, x_T, x_F);
fhhos4::JoinVector(x_T, x_F, cells, faces, x); // the DoFs of x in no list unchanged (e.g. the Dirichlet values)
```

- `M`: the whole matrix (both lower and upper triangular parts), square, an Eigen row-major matrix (or a map of your
  arrays) with the index type `std::int32_t` or `std::int64_t`.
- The lists: `std::int32_t` or `std::int64_t`; each DoF in one list at most. The order of the lists gives the order of
  the blocks, and that of the DoFs of each cell or face the order of its basis (reorder them, e.g. by degree, here).
- The Dirichlet DoFs, in no list: move their contribution to the right-hand side first (`b - M x_D`).
- `ExtractBlocksSym()` gives the lower triangular parts of `A_TT` and `A_FF` (what `SetupSym()` reads),
  `ExtractBlocksNonSym()` the four blocks in full.

[With MPI](#with-mpi), from blocks stored in two parts (e.g. Eigen maps of the CSR arrays of PETSc's
`MatMPIAIJGetSeqAIJ()`, or of hypre's `diag` and `offd`):

```cpp
// diag: the columns of the faces of the process; offd: the others, compressed, whose global DoF numbers are
// offdColumns (PETSc's garray, hypre's col_map_offd)
fhhos4::DistributedBlocks b = fhhos4::MergeOffDiag(A_TF_diag, A_TF_offd, A_TF_offdColumns, A_FF_diag, A_FF_offd, A_FF_offdColumns, faceBlockSize);
solver.SetupSym(A_TT_diag, b.A_TF, b.A_FF, cellInterpOfOne, faceInterpOfOne, b.FaceGhostGlobalNumbers);
```

- The global DoF numbers: the DoF k of the face f is `f * faceBlockSize + k`, the faces numbered in rank order (as
  `faceGhostGlobalNumbers`): those of PETSc's and hypre's matrices when each process owns the DoFs of whole faces.
  `offdColumns`: strictly increasing.
- The parts: Eigen sparse matrices, compressed, row-major, all with the same index type `std::int32_t` or
  `std::int64_t`; `offdColumns`: pointers to `std::int32_t` or `std::int64_t`.
- The local columns: the faces of the process, then its ghosts (those of the off-diagonal columns of `A_TF` and
  `A_FF`) by increasing global number, in `b.FaceGhostGlobalNumbers`.
- `A_TT` has no off-diagonal part (cells are not shared): give its diagonal part to the setup, read in place if its
  index type is the library's.
- For `SetupNonSym()`, the overload with `A_FT`: `MergeOffDiag(A_TF_diag, ..., A_FT_diag, A_FT_offd,
  A_FT_offdColumns, A_FF_diag, ..., faceBlockSize, cellBlockSize)`. The off-diagonal columns of `A_FT` are the cells
  of other processes (the DoF k of the cell c: `c * cellBlockSize + k`). It gives `b.A_FT` and
  `b.CellGhostGlobalNumbers` too.

### Parameters

The public members of `fhhos4::Solver`, read by the setup (`Tolerance` and `MaxIterations` by `Solve()`), with the
library's defaults (`library/include/fhhos4.h`, which details each one):

| Parameter | Default | Values |
|---|---|---|
| `PolynomialSpace` | `PolynomialSpace::P` | Space of the cell and face bases: <ul><li>`P`: degree at most k</li><li>`Q`: degree at most k in each variable (quadrilaterals, hexahedra)</li></ul> |
| `VectorComponents` | `1` | `>= 1`: components of the unknown, interleaved in each block (see [Inputs](#inputs)) |
| `SingularMatrix` | `false` | `true`: the kernel is the interpolation of 1 of each component (pure Neumann; translations, not rotations). `SetupSym()` only. The solution is defined up to that kernel |
| `Krylov` | `FCG` if symmetric, `GCR` if non-symmetric | <ul><li>`FCG`: flexible conjugate gradient</li><li>`GCR`: restarted GCR, any matrix</li><li>`None`: multigrid cycles alone (e.g. [as a preconditioner](#use-as-a-preconditioner))</li></ul> |
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
| `CoarseMatrixMaxRows` | `1000` | Maximum number of rows of the coarsest matrix (on several processes: in total) |
| `MinRowsPerProcess` | (MPI only) `50` | On several processes, the coarsening also stops at this many rows per process on average (`0`: only `CoarseMatrixMaxRows`). Each process aggregates its own cells: on many processes, the coarsest matrix cannot get down to `CoarseMatrixMaxRows` |
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
