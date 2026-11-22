[← The fhhos4 library](../../README_library.md)

# The system to give

## Inputs

- **The blocks** $A_{TT}$, $A_{TF}$, $A_{FF}$, and $A_{FT}$ for
  a non-symmetric matrix. The DoFs of each cell (resp. face) are contiguous, cells and faces in any order. In a
  symmetric system on one process, $A_{TT}$ and $A_{FF}$ are read from their lower triangular part (lower-only or full);
  otherwise in full ($A_{FF}$ on several processes, all the blocks of a non-symmetric system). A triangular block where
  both lower and upper triangular parts are read (or an upper triangular one where the lower is) is refused.
- **The format of the matrices**: Eigen sparse matrices, compressed, row-major, with the index type of the library
  (`fhhos4::SparseIndex`: `int`, or 64 bits with `-DENABLE_64BIT_INDEX=ON`): `fhhos4::SparseMatrix`, or an Eigen map of your
  own CSR arrays, `fhhos4::SparseMatrixMap(rows, cols, nonZeros, rowPointers, columnIndices, values)`, the columns
  strictly increasing in each row (sorted, no duplicates; checked). **The matrices are read in place, never copied nor
  converted**: assemble them in this format. Any other Eigen sparse matrix converts by Eigen, in one line, a copy that
  stays visible in your code: `fhhos4::SparseMatrix B = A;` (see below).
- **The near null space** of the operator, interpolated on the cell and face bases: `cellNearNullSpace` and
  `faceNearNullSpace`, one row per DoF (the rows of $A_{TF}$ for the cells, its columns for the faces) and one column per
  function, the same functions in both. It is what your own interpolation gives, and what PETSc's
  `MatSetNearNullSpace()` takes. An Eigen vector is a matrix of one column.
  - **A scalar unknown**: one column, the interpolation of the function 1.
  - **A vector-valued unknown** of $n$ components: $n$ columns, the interpolations of the constants $e_1, \dots, e_n$,
    in this order.
  - **Linear elasticity**: the rigid body modes, i.e. the translations and the rotations (3 columns in 2D, 6 in 3D),
    in any order, or any basis of the space they span. Without the rotations, the solver only converges if the
    displacement is imposed on the whole boundary. With them: `SetupSym()`, `FaceDegree` and `CellDegree` at least
    1; with `MixedPrecision`, only the levels of degree >= 1 are in single precision. The pure traction problem (no
    Dirichlet face): `SingularMatrix`, whose kernel is then the rigid body modes given.

  The bases must be **hierarchical** (ordered by degree: for each d, the first functions of a block span P^d, or Q^d
  whose degree is the largest partial degree; with `PCoarsening::ToZero`, only the degree 0 first, or the degrees 0 and
  1 for the rotations) **with a constant first function** $\varphi_0$, of any scale, which may vary from one cell or face
  to another (e.g. L2-orthonormal bases). The setup checks it on the near null space: the interpolation of a constant
  has one non-zero coefficient per cell and per face, the first (of its component).
- **The dimension and the degrees**: set `Dimension` (2 or 3), `FaceDegree` (k), `CellDegree` (l: k - 1, k or k + 1,
  refused otherwise) and `PolynomialSpace` (P or Q). The setup checks the sizes and structure of the inputs against them
  (full polynomial spaces: P^l in the cells and P^k on the faces, or Q^l and Q^k; times `VectorComponents`).
- **Left unset** (`-1`, `PolynomialSpace::Inferred`: the defaults), they are inferred by the setup from the blocks,
  which is handy to try the library: from the stored entries of $A_{TT}$ and $A_{TF}$, the setup reads the DoFs per cell and
  per face, and finds the discretizations that have them. **Set `VectorComponents` then** (`1` for a scalar unknown):
  the DoFs don't tell it. It assumes
  - the same degrees in every cell and face, `CellDegree` k - 1, k or k + 1;
  - **the zeros inside the blocks of each cell stored**: the setup reads where the entries are, not their values;
  - a face between two cells (on several processes: of the same process).

  If different discretizations have these DoFs and would need different solvers (e.g. 64 DoFs per cell and 9 per
  face: Q in 3D with k = 2 and l = 3, or in 2D with k = 8 and l = 7), the setup stops with an error that names the
  parameter to set; so it does if none has them. Otherwise, it writes the values found into these parameters (printed
  with `Verbosity = 1`). With `Dimension` and the degrees set, `PolynomialSpace` alone is inferred, from the stored
  entries too (P if they don't tell), and `VectorComponents` is 1 if not set.
- **A vector-valued unknown** of $n$ components: `VectorComponents = n`. Each cell and face block holds the scalar basis
  times each unit vector, **interleaved by component** (array of structures): $\varphi_i e_c$ is the DoF $n i + c$ of the
  block.

## Other matrix formats

`fhhos4` provides the following helper functions to convert your matrices to the format `fhhos4` expects (see
[Inputs](#inputs)).

| Your matrices | How to convert | Makes a copy? |
|---|---|---|
| Column-major (Eigen's default), or another index type | `fhhos4::SparseMatrix B = A;` (Eigen's conversion) | Yes |
| Column-major and symmetric ($A_{TT}$, $A_{FF}$, or $A$ of a symmetric system) | `TransposedView(A)`: $A^T$ as an Eigen map of the arrays of `A`, which is $A$, up to the rounding errors of its assembly | No |
| Your own CSR arrays, with the index type of the library | `fhhos4::SparseMatrixMap(rows, cols, nonZeros, rowPointers, columnIndices, values)` | No |
| The monolithic matrix of the system (the cells and faces DoFs in your numbering, the Dirichlet DoFs included) | `ExtractBlocksSym(M, cells, faces)`, `ExtractBlocksNonSym(M, cells, faces)`; `SplitVector()` and `JoinVector()` for the vectors | Yes |
| With MPI, blocks stored in two parts, as PETSc (MPIAIJ) and hypre (ParCSR) store them: the columns of the faces of the process, and the others with their global numbers | `MergeOffDiag(...)` (below) | Yes (not $A_{TT}$) |

The copies are Eigen matrices of your code (`fhhos4::SparseMatrix`). From a monolithic matrix (e.g. assembled by a
framework in its own numbering):

```cpp
// The DoFs of each cell (resp. face) in the numbering of M: those of the i-th are DoFs[Offsets[i]], ...,
// DoFs[Offsets[i + 1] - 1], in the order of its basis; the DoFs in no list (e.g. of the Dirichlet faces) are dropped
fhhos4::DoFLists<std::int64_t> cells{ nCells, cellOffsets, cellDoFs }, faces{ nFaces, faceOffsets, faceDoFs };
Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor, std::int64_t>> M(n, n, nnz, rowPointers, columns, values);
fhhos4::HybridBlocks blocks = fhhos4::ExtractBlocksSym(M, cells, faces); // A_TT, A_TF, A_FF: one parallel pass
solver.SetupSym(blocks.A_TT, blocks.A_TF, blocks.A_FF, cellNearNullSpace, faceNearNullSpace);
fhhos4::SplitVector(b, cells, faces, b_T, b_F); // b: in the numbering of M
solver.Solve(b_T, b_F, x_T, x_F);
fhhos4::JoinVector(x_T, x_F, cells, faces, x); // the DoFs of x in no list unchanged (e.g. the Dirichlet values)
```

- `M`: the whole matrix (both lower and upper triangular parts), square, an Eigen row-major matrix (or a map of your
  arrays) with the index type `std::int32_t` or `std::int64_t`.
- The lists: `std::int32_t` or `std::int64_t`; each DoF in one list at most. The order of the lists gives the order of
  the blocks, and that of the DoFs of each cell or face the order of its basis (reorder them, e.g. by degree, here).
- The Dirichlet DoFs, in no list: move their contribution to the right-hand side first ($b - M x_D$).
- `ExtractBlocksSym()` gives the lower triangular parts of $A_{TT}$ and $A_{FF}$ (what `SetupSym()` reads),
  `ExtractBlocksNonSym()` the four blocks in full.

[With MPI](mpi.md), from blocks stored in two parts (e.g. Eigen maps of the CSR arrays of PETSc's
`MatMPIAIJGetSeqAIJ()`, or of hypre's `diag` and `offd`):

```cpp
// diag: the columns of the faces of the process; offd: the others, compressed, whose global DoF numbers are
// offdColumns (PETSc's garray, hypre's col_map_offd)
fhhos4::DistributedBlocks b = fhhos4::MergeOffDiag(A_TF_diag, A_TF_offd, A_TF_offdColumns, A_FF_diag, A_FF_offd, A_FF_offdColumns, faceBlockSize);
solver.SetupSym(A_TT_diag, b.A_TF, b.A_FF, cellNearNullSpace, faceNearNullSpace, b.FaceGhostGlobalNumbers);
```

- The global DoF numbers: the DoF k of the face f is `f * faceBlockSize + k`, the faces numbered in rank order (as
  `faceGhostGlobalNumbers`): those of PETSc's and hypre's matrices when each process owns the DoFs of whole faces.
  `offdColumns`: strictly increasing.
- The parts: Eigen sparse matrices, compressed, row-major, all with the same index type `std::int32_t` or
  `std::int64_t`; `offdColumns`: pointers to `std::int32_t` or `std::int64_t`.
- The local columns: the faces of the process, then its ghosts (those of the off-diagonal columns of $A_{TF}$ and
  $A_{FF}$) by increasing global number, in `b.FaceGhostGlobalNumbers`.
- $A_{TT}$ has no off-diagonal part (cells are not shared): give its diagonal part to the setup, read in place if its
  index type is the library's.
- For `SetupNonSym()`, the overload with $A_{FT}$: `MergeOffDiag(A_TF_diag, ..., A_FT_diag, A_FT_offd,
  A_FT_offdColumns, A_FF_diag, ..., faceBlockSize, cellBlockSize)`. The off-diagonal columns of $A_{FT}$ are the cells
  of other processes (the DoF k of the cell c: `c * cellBlockSize + k`). It gives `b.A_FT` and
  `b.CellGhostGlobalNumbers` too.
