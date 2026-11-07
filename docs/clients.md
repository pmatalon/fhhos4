# The library's clients and targets

How each client code (or target) assembles its matrix, whether it runs on several processes, and how it reaches the
library: what drives the library's interface and helpers. Keep it up to date when an integration changes or a target
is studied. Terms: `GLOSSARY.md` (blocks, monolithic matrix, distributed matrix). Facts about external codes are
from their sources or documentation unless marked "unchecked".

| Code | Status | Who assembles, in what format | What it gives | MPI | How it reaches the library |
|---|---|---|---|---|---|
| The program | in the repository (`-s libuamg`) | itself, Eigen CSR (triplets) | blocks | yes (GMSH parts) | `SetupSym()`/`SetupNonSym()`; on several processes `DistributedHybridSystem` copies into the library's layout (no copy planned) |
| HArDCore3D | client (fork, branch `fhhos4`) | itself, Eigen row-major | blocks | no | `SetupSym()` (diffusion), `SetupNonSym()` (diffusion-advection-reaction) |
| Disk++ | client (fork, branch `fhhos4`) | itself, Eigen row-major | blocks | no | `SetupSym()`, `VectorComponents` (vector Laplacian) |
| NGSolve | client (the adapter, Python) | itself (`ngla`), CSR; partial rows on several processes | monolithic matrix (product space) | optional | `ExtractBlocksSym()`/`NonSym()`, `SplitVector()`/`JoinVector()`; on several processes the adapter sums the shared rows itself |
| PETSc users | target: the preconditioner (`mpi/12`) | PETSc (`MatSetValues`, `MatAssemblyEnd`) | MATNEST of AIJ blocks, or monolithic AIJ + index sets | yes (or one process) | merge helper (to write); `MatCreateSubMatrix()` then the merge helper for a monolithic matrix; `SplitVector()`/`JoinVector()` |
| DOLFINx | target (`mpi/04`) | itself (`MatrixCSR`) or PETSc (AIJ, MATNEST) | blocks or monolithic matrix | always | through PETSc: as PETSc users; native `MatrixCSR`: to study |
| Code_Aster | target, studied only (no code yet) | itself (`MATR_ASSE`, MORSE storage); PETSc in its parallel versions | monolithic matrix, uncondensed | yes | one process: `ExtractBlocksSym()`/`NonSym()`, `SplitVector()`/`JoinVector()`; several: through its PETSc conversion (to study) |
| MFEM | possible | itself, then hypre (`ParallelAssemble()`: `HypreParMatrix`, ParCSR) | blocks (`BlockOperator`) or monolithic matrix | yes | merge helper (ParCSR has PETSc's layout); unchecked |
| deal.II | possible | itself, into PETSc or Trilinos (Tpetra) matrices | blocks (block matrices) or monolithic matrix | yes | through PETSc: as PETSc users; Tpetra (one CSR, local columns, any global numbers): no helper; unchecked |

## Notes by code

**The program.** `Diffusion_HHO` assembles `A_T_T`, `A_T_ndF` and `A_ndF_ndF` (in full), the Dirichlet faces apart.
On several processes it recomputes a layer of ghost cells, then copies its blocks into the library's layout
(`SolveDistributed()`); planned (ADR-0011, `.scratch/mpi/issues/12-petsc-binding.md`): its own cells only, its own
cells and faces first in its local numbering, the blocks given as their first rows (no copy).

**HArDCore3D.** Assembles the four blocks and the right-hand sides in its own loops, the columns of the Dirichlet faces
apart (moved to the right-hand sides). L2-orthonormal bases: the interpolations of 1 are not ones.

**Disk++.** `assemble_uncondensed()` beside its condensed assembly (`poisson_solver`, `vector_laplacian`). Hierarchical
monomials, phi_0 = 1. Its condensed matrix is symmetric only to rounding.

**NGSolve.** HDG is one bilinear form on a product space: a monolithic matrix (uint64 row pointers, int32 columns, read
in place), its DoFs sorted by degree by the adapter (Dubiner order). On several processes `local_mat` holds partial
rows (each process's own cells only, no global numbering): the adapter extracts the blocks per process, numbers the
faces (exclusive scan) and sums the rows of the shared faces on their owners (`.scratch/ngsolve/issues/03-...`).
Non-symmetric: one process only. Towards hypre and PETSc NGSolve sends its partial rows and lets them sum
(`HYPRE_IJMatrixAddToValues`; ngsPETSc: MATIS, then AIJ).

**PETSc users.** Each process assembles its own cells; PETSc sends the contributions to rows of other processes. A
block on several processes is an MPIAIJ: two CSRs (the owned and the other columns) and `garray`, the global DoF
numbers of the other columns. A copy before the library is accepted (as for hypre, AmgX, MUMPS). The interpolation of
1 is the near-null space vector (`MatGetNearNullSpace`).

**DOLFINx.** Each process assembles its own cells into owned and ghost rows (local columns: owned, then ghosts), then
`scatter_rev()` adds the ghost rows to their owners; or PETSc matrices, including MATNEST (nest assembly). HDG: a cell
space and a facet space (mixed dimensions); its demos condense in the kernel, which the library cannot take
(ADR-0002): the uncondensed blocks must be assembled. Whether the column index maps of `A_TF` and `A_FF` coincide in a
native `MatrixCSR`: unchecked.

**Code_Aster** (documentation, 2026-11-07). HHO: mechanics `3D_HHO`, `D_PLAN_HHO` (and `_GRAD_HHO`), non-linear
operators only (`STAT_NON_LINE`, `DYNA_NON_LINE`); thermal `PLAN_HHO`, `AXIS_HHO`, `3D_HHO`, degrees 0 to 4. Its HHO
has no static condensation ("la condensation statique n'est pas implémentée", R3.06.14): it assembles the uncondensed
system, cell and face unknowns, into one matrix. Mechanics: a vector-valued unknown (`VectorComponents`), a new matrix
at each Newton iteration (a setup each time). Several processes: the elements distributed over the processes
(`MATR_DISTRIBUEE`), PETSc only in the parallel versions, with its own global numbering (blocks of contiguous rows per
process); its preconditioners through PETSc (`PRE_COND`: `LDLT_SP`, `GAMG`, `BOOMER`, `ML`, `HPDDM`, `FIELDSPLIT`,
and `UTILISATEUR`, a preconditioner in Python: a possible entry point, unchecked). Unchecked: whether MORSE stores half
of a symmetric matrix; whether its Dirichlet conditions by Lagrange multipliers (`AFFE_CHAR_MECA`) can be avoided by
elimination (`AFFE_CHAR_CINE`) with HHO: the multigrid preconditioners require no multipliers.

**MFEM, deal.II.** Not studied: from memory.

## What it means for the library

- Blocks assembled by the client (the program, HArDCore3D, Disk++): the setups alone, no helper.
- A monolithic matrix is common (NGSolve, Code_Aster, PETSc's AIJ, DOLFINx's monolithic assembly): `ExtractBlocks*()`
  on one process; on several processes, through PETSc (`MatCreateSubMatrix()`) or, for a client without PETSc, a
  distributed `ExtractBlocks*()` (not planned).
- PETSc's and hypre's distributed matrices (two CSRs per block): the merge helper (`.scratch/mpi/issues/12-...`).
- Partial rows (NGSolve, PETSc's MATIS): summed by the client today.
- Repeated setups on the same pattern (Code_Aster's Newton iterations): the setup time matters.
