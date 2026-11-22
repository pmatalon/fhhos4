[← The fhhos4 library](../../README_library.md)

# Python

The Python module `fhhos4` has the classes, functions, parameters and enums of the [C++ library](cpp.md), with the same names
(`solver.Krylov = fhhos4.Krylov.GCR`), and also reads its inputs in place, as NumPy arrays.

- [Installation](#installation)
- [How to use](#how-to-use)
  - [What differs from C++](#what-differs-from-c)
- Common to the three interfaces: [The system to give](system.md), [Several processes (MPI)](mpi.md),
  [Parameters](parameters.md)

## Installation

It needs Python (>= 3.9),
[nanobind](https://nanobind.readthedocs.io) and NumPy (SciPy is optional), and is built with the library (not installed
yet: use it from the build tree):

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_PYTHON=ON
cmake --build build
export PYTHONPATH=$PWD/build/python    # the module fhhos4
ctest --test-dir build                 # its tests (they need SciPy)
```

## How to use

```python
import numpy as np
import fhhos4

solver = fhhos4.Solver()
solver.VectorComponents = 1  # 1 if scalar; the next four optional: inferred from the blocks if not set
solver.Dimension = 3
solver.FaceDegree = k
solver.CellDegree = k
solver.PolynomialSpace = fhhos4.PolynomialSpace.P
# e.g. scipy.sparse.csr_array; the near null space: NumPy vectors, or arrays of one column per function in Fortran order
solver.SetupSym(A_TT, A_TF, A_FF, cellNearNullSpace, faceNearNullSpace)
x_T, x_F = np.zeros(b_T.size), np.zeros(b_F.size)
result = solver.Solve(b_T, b_F, x_T, x_F)  # x_T, x_F written in place
```

### What differs from C++

- **The matrices**: SciPy's CSR matrices (or any object with the attributes `indptr`, `indices`, `data`, `shape` and
  `format` `"csr"`), a `fhhos4.CSRView(rows, cols, rowPointers, columnIndices, values)` of NumPy arrays (or of any
  contiguous array with the buffer protocol), or a `fhhos4.CSRMatrix` (what the helpers return). The values: float64;
  the columns sorted (SciPy: `sort_indices()`, in place). The setups take the index type of the library (int32 by
  default) for both index arrays, the helpers int32 or int64 (uint64 read as int64). Other formats raise `TypeError`:
  convert them with `fhhos4.ToCSR(M)` (a copy, SciPy's CSC matrices included), or view a symmetric CSC matrix with
  `fhhos4.TransposedView(A)`.
- **The vectors**: NumPy vectors of float64, read (outputs: written) in place; another dtype raises `TypeError`. The
  lists of `fhhos4.DoFLists(count, offsets, dofs)`: NumPy arrays of the same integer dtype.
- **Lifetimes**: the solver keeps alive the arrays that its solves read (those of $A_{FT}$, for `SetupNonSym()`).
  `CondensedMatrix()` returns a copy. The arrays of a `fhhos4.CSRMatrix` (`RowPointers`, `ColumnIndices`, `Values`)
  are read-only NumPy views of the library's memory, which they keep alive.
- **The errors** raise `fhhos4.Error`, a `RuntimeError`. The setups, solves and helpers release the GIL.
- **With MPI** (`-DENABLE_MPI=ON` too; `fhhos4.WithMPI` is then `True`): `solver.Comm` is an
  [mpi4py](https://mpi4py.readthedocs.io) communicator (`MPI.COMM_WORLD`, or its Fortran handle `comm.py2f()`). Import
  mpi4py, built with the same MPI as the library, before `fhhos4` (it initializes MPI). The ghost numbers are NumPy
  arrays of int64, read in place (`None` by default: no ghosts). Without `Comm`, mpi4py is not needed. The
  `offdColumns` of `MergeOffDiag()`: NumPy arrays of int32 or int64; its ghosts: arrays of int64.
