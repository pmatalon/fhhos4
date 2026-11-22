[← The fhhos4 library](../../README_library.md)

# PETSc

`fhhos4` is also a preconditioner of [PETSc](https://petsc.org), `-pc_type fhhos4`, for its Krylov methods on the
uncondensed system.

- [Installation](#installation)
  - [Loading it](#loading-it)
- [How to use](#how-to-use)
  - [The operator](#the-operator)
  - [The discretization](#the-discretization)
  - [The near null space](#the-near-null-space)
  - [Symmetric or not](#symmetric-or-not)
  - [The Krylov method](#the-krylov-method)
  - [Options](#options)
- Common to the three interfaces: [The system to give](system.md), [Several processes (MPI)](mpi.md),
  [Parameters](parameters.md)

## Installation

Build it with `-DENABLE_PETSC_PC=ON`, against the PETSc of your code (real scalars in double
precision; with MPI, `-DENABLE_MPI=ON` and the same MPI).

### Loading it

Without recompiling, `-dll_append <prefix>/lib/libfhhos4_petsc.so` among PETSc's options (command
line, `PETSC_OPTIONS`, `petsc4py.init()`); or link `fhhos4::petsc` and call `PetscDLLibraryRegister_fhhos4_petsc()`
once after `PetscInitialize()` (`#include <fhhos4_petsc.h>`, which defines `PCFHHOS4`).

## How to use

```bash
mpiexec -n 4 ./my_code -dll_append <prefix>/lib/libfhhos4_petsc.so -ksp_type fcg -pc_type fhhos4 \
    -pc_fhhos4_vector_components 1 -pc_fhhos4_dimension 3 -pc_fhhos4_face_degree 1 -pc_fhhos4_cell_degree 1
```

Complete examples, an HHO discretization of degree 1 or 2 assembled cell by cell, its interpolation of 1, and the
solve: [`examples/petsc_preconditioner.cpp`](../../examples/petsc_preconditioner.cpp) (C++, linked to `fhhos4::petsc`), and
[`examples/petsc_preconditioner.py`](../../examples/petsc_preconditioner.py) (petsc4py, loading the library by `-dll_append`).

### The operator

A `MATNEST` of 2 x 2 AIJ blocks `[[A_TT, A_TF], [A_FT, A_FF]]`, the cells first, assembled as usual
(`MatSetValues()` of the local matrices of each cell, `MatAssemblyEnd()`); or one AIJ matrix in your numbering, with
the rows of its cells and of its faces, `PCFHHOS4SetIS(pc, cells, faces)` before the setup (petsc4py:
`pc.compose("fhhos4_cells", cells)` and `pc.compose("fhhos4_faces", faces)`), listed cell by cell and face by face,
the DoFs of each in the order of [Inputs](system.md#inputs) (the components interleaved). The rows in neither index set (e.g.
Dirichlet DoFs eliminated by `MatZeroRowsColumns()`) get the identity. On several processes, each process owns the
DoFs of whole cells and faces (and gives its own rows in the index sets). The setup copies the blocks into the format
of the library (`MergeOffDiag()`, see [Other matrix formats](system.md#other-matrix-formats)) but $A_{TT}$, read in place; from
one AIJ matrix, it first extracts them (`MatCreateSubMatrix()`, freed after the setup). The vectors: plain
(`MatCreateVecs()` of the operator) or, with a MATNEST, nested (`VECNEST`, the cells then the faces).

### The discretization

Give its options (the table below); or, to try the preconditioner, only
`-pc_fhhos4_vector_components` (`1` for a scalar unknown): the others are then inferred from the blocks, under the
assumptions of [Inputs](system.md#inputs), and an ambiguity names the option to give. With the dimension and the degrees
given, the space too is read from the blocks, unless P and Q have the same DoFs (e.g. at degree 0); P if the blocks
don't tell.

### The near null space

The near null space of the operator (no default), `MatSetNearNullSpace()`: the vectors of [Inputs](system.md#inputs), on
all the DoFs (a scalar unknown: the interpolation of 1 on your bases; a vector-valued one: of the constant of each
component; elasticity: the rigid body modes interpolated on your bases, orthonormalized as PETSc requires: any
basis of their span is accepted. Not `MatNullSpaceCreateRigidBody()`, whose vectors are nodal values).
At degree 0 on one component, it can be the constant vector (`MatNullSpaceCreate(comm, PETSC_TRUE, 0, NULL, ...)`).

### Symmetric or not

The symmetric setup if the operator is known to be symmetric
(`MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE)`, or `MAT_SPD`) or with `-pc_fhhos4_symmetric`, the non-symmetric one
otherwise.

### The Krylov method

Each application is one K-cycle, a non-linear preconditioner: use a flexible method, `-ksp_type
fcg` (symmetric), `fgmres` or `gcr`, not PETSc's default GMRES. `-ksp_view` shows the setup. PETSc's FCG tests
the preconditioned residual by default, far from the true one on badly conditioned systems (elasticity of thin
structures): `-ksp_norm_type unpreconditioned` for `-ksp_rtol` on $b - Ax$, as the library's own solve.

### Options

| Option | Default | Parameter |
|---|---|---|
| `-pc_fhhos4_vector_components` | none: required if the dimension or a degree is not given; `1` otherwise | `VectorComponents` |
| `-pc_fhhos4_dimension`, `-pc_fhhos4_face_degree`, `-pc_fhhos4_cell_degree` | inferred from the blocks (see [Inputs](system.md#inputs)) | `Dimension`, `FaceDegree`, `CellDegree` |
| `-pc_fhhos4_polynomial_space p\|q` | inferred from the blocks | `PolynomialSpace` |
| `-pc_fhhos4_symmetric` | from the operator | the symmetric setup |
| `-pc_fhhos4_cycle k\|v\|w` | `k` | `Cycle` |
| `-pc_fhhos4_krylov none\|fcg\|gcr` | `none`: one cycle | `Krylov` in each application: the whole solve of the library, with `-ksp_type preonly` |
| `-pc_fhhos4_max_it`, `-pc_fhhos4_rtol` | `1`, `1e-8` | `MaxIterations` (cycles with `none`), `Tolerance` |
| `-pc_fhhos4_verbosity` | `0` | `Verbosity` |

The other parameters keep the library's defaults.
