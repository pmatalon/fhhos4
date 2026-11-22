[← The fhhos4 library](../../README_library.md)

# Several processes (MPI)

If the library is compiled with MPI (`-DENABLE_MPI=ON`), `solver.Comm` sets the processes that solve the system:

- **Not set** (`MPI_COMM_NULL`, the default) or `MPI_COMM_SELF`: each process solves its own system alone, as without
  MPI. If `Comm` is not set, the solver makes no MPI call: a code that doesn't use MPI needs no change, and doesn't even
  have to initialize MPI (unless it uses `CoarseSolver::MUMPS`).
- **A communicator of several processes**: the processes solve one system together, distributed by rows. Every setup
  and solve is collective: all the processes of the communicator call it, with the same parameters.

**Who owns what.** Each cell and each face belongs to one process. A face must belong to a process that owns one of
its cells (a face has one or two cells), e.g. the one of smaller rank. A process may own no cell and no face: it then
gives empty matrices and vectors, and still calls the setups and solves.

## Expected format

Each process gives `SetupSym()` the following (the `FHHOS4_MPI` branch of the example above):

- $A_{TT}$ and $A_{TF}$: the rows of its cells. $A_{FF}$: the rows of its faces, assembled (with the contributions of the
  cells of other processes) and with both lower and upper triangular parts (on one process, the lower triangular part
  alone is enough).
- The columns of $A_{TF}$ and $A_{FF}$ follow a local numbering of the faces: first the faces of the process, in the order
  of the rows of $A_{FF}$, then its **ghost faces**, in any order. The ghost faces are the faces of other processes that
  appear in these rows: the faces of its cells that other processes own, and the faces coupled to its own faces in
  $A_{FF}$ through the cells of other processes. This is the local numbering of DOLFINx's index maps.
- `faceGhostGlobalNumbers`: the global numbers of its ghost faces (faces, not DoFs), in the order of their local
  numbers. The faces are numbered as PETSc and hypre require for their distributed matrices: those of each process
  contiguously, in rank order.
- The vectors $b_T$, $b_F$, $x_T$, $x_F$ and the near null space: the values on its own cells and faces.

For example, with two processes, rank 0 owns the faces 0 to 9 and rank 1 the faces 10 to 19. If the rows of rank 0
contain the faces 15 and 12 of rank 1, rank 0 numbers them locally 10 and 11 (or 11 and 10) and gives
`faceGhostGlobalNumbers = {15, 12}` (or `{12, 15}`).

[`examples/fhhos4_solver.cpp`](../../examples/fhhos4_solver.cpp) assembles the rows of each process in this format, its
ghost faces numbered in the order it meets them.

PETSc's and hypre's distributed matrices store each block in two parts: `MergeOffDiag()` converts them (see
[Other matrix formats](system.md#other-matrix-formats)).

## Non-symmetric systems

Each process gives the rows of its faces in $A_{FT}$. Its columns follow a local numbering of
the cells, as those of $A_{TF}$ and $A_{FF}$ for the faces: first the cells of the process, in the order of the rows of
$A_{TT}$, then its **ghost cells**, in any order (the cells of other processes that appear in these rows). $A_{FT}$ is
the only block with ghost cells. Their global numbers, `cellGhostGlobalNumbers`, number the cells as
`faceGhostGlobalNumbers` the faces:
```cpp
solver.SetupNonSym(A_TT, A_TF, A_FT, A_FF, cellNearNullSpace, faceNearNullSpace, faceGhostGlobalNumbers,
                   cellGhostGlobalNumbers);
```
On each process, the Gauss-Seidel smoothers follow the downwind order of its own faces.

## What you get back

`Solve()` writes $x_T$ and $x_F$ on the cells and faces of the process, and returns the `Result`
of the whole system. An error in the inputs of one process is raised on all the processes, so that none is left
waiting.

## What changes in the solver

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
  (its matrix is small). The same holds for `CoarseSolver::MUMPS` on one process. With the OpenMP build of OpenBLAS,
  if your code sets its threads (`openblas_set_num_threads()`, which the environment can't override), its products
  inside MUMPS hang: call `openblas_set_num_threads(1)` before the setups and solves.

In Python: see [Python](python.md).
