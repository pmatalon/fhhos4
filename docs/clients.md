# The library's clients and targets

How each client code (or target) assembles its matrix, whether it runs on several processes, and how it reaches the
library: what drives the library's interface and helpers. Keep it up to date when an integration changes or a target
is studied. Terms: `GLOSSARY.md` (blocks, monolithic matrix, distributed matrix). Facts about external codes are
from their sources or documentation unless marked "unchecked".

| Code | Status | Who assembles, in what format | What it gives | MPI | How it reaches the library |
|---|---|---|---|---|---|
| The program | in the repository (`-s libuamg`) | itself, Eigen CSR (triplets) | blocks | yes (GMSH parts) | `SetupSym()`/`SetupNonSym()`; on several processes `DistributedHybridSystem` copies into the library's layout (no copy planned) |
| HArDCore3D | client (fork, branch `fhhos4`) | itself, Eigen row-major | blocks | no | `SetupSym()` (diffusion), `SetupNonSym()` (diffusion-advection-reaction) |
| Disk++ | client (fork, branch `fhhos4`) | itself, Eigen row-major | blocks | no | `SetupSym()`, `VectorComponents` (vector Laplacian), the rigid body modes (linear elasticity) |
| NGSolve | client (the adapter, Python) | itself (`ngla`), CSR; partial rows on several processes | monolithic matrix (product space) | optional | `ExtractBlocksSym()`/`NonSym()`, `SplitVector()`/`JoinVector()`; on several processes the adapter sums the shared rows itself |
| PETSc users | target: the preconditioner (`mpi/12`) | PETSc (`MatSetValues`, `MatAssemblyEnd`) | MATNEST of AIJ blocks, or monolithic AIJ + index sets | yes (or one process) | `-pc_type fhhos4` (`ENABLE_PETSC_PC`: `MergeOffDiag()` from `MatMPIAIJGetSeqAIJ()`); `MatCreateSubMatrix()` then the merge for a monolithic matrix; `SplitVector()`/`JoinVector()` |
| DOLFINx | target (`mpi/04`) | itself (`MatrixCSR`) or PETSc (AIJ, MATNEST) | blocks or monolithic matrix | always | through PETSc: as PETSc users; native `MatrixCSR`: to study |
| Code_Aster | target, studied only (no code yet) | itself (`MATR_ASSE`, MORSE storage); PETSc in its parallel versions | monolithic matrix, uncondensed | yes | through PETSc: a `KSP_UTIL` function of the user's script (petsc4py) builds the index sets, the blocks and `-pc_type fhhos4` (see below) |
| MFEM | possible | itself, then hypre (`ParallelAssemble()`: `HypreParMatrix`, ParCSR) | blocks (`BlockOperator`) or monolithic matrix | yes | `MergeOffDiag()` (ParCSR has PETSc's layout); unchecked |
| deal.II | possible | itself, into PETSc or Trilinos (Tpetra) matrices | blocks (block matrices) or monolithic matrix | yes | through PETSc: as PETSc users; Tpetra (one CSR, local columns, any global numbers): no helper; unchecked |

## Notes by code

**The program.** `Diffusion_HHO` assembles `A_T_T`, `A_T_ndF` and `A_ndF_ndF` (in full), the Dirichlet faces apart.
On several processes it recomputes a layer of ghost cells, then copies its blocks into the library's layout
(`SolveDistributed()`); planned (ADR-0011, `.scratch/mpi/issues/12-petsc-binding.md`): its own cells only, its own
cells and faces first in its local numbering, the blocks given as their first rows (no copy).

**HArDCore3D.** Assembles the four blocks and the right-hand sides in its own loops, the columns of the Dirichlet faces
apart (moved to the right-hand sides). L2-orthonormal bases: the interpolations of 1 are not ones.

**Disk++.** `assemble_uncondensed()` beside its condensed assembly (`poisson_solver`, `vector_laplacian`,
`linear_elasticity`). Hierarchical monomials, phi_0 = 1. Its condensed matrix is symmetric only to rounding.
`linear_elasticity` (Di Pietro-Ern, k >= 1) gives the near null space by its own projection (`project_function()`) of
the translations and, with `-r`, of the rotations about the origin, on each cell and each face in the system; its
vectors pass the library's checks. `-t`: pure traction (the traction of the manufactured solution on the whole
boundary of the unit square or cube, no Dirichlet face): `SingularMatrix` with the rigid body modes, the L2 error after
the removal of the rigid motion (its L2 projection on the modes). The only elasticity the library is run on.

**NGSolve.** HDG is one bilinear form on a product space: a monolithic matrix (uint64 row pointers, int32 columns, read
in place), its DoFs sorted by degree by the adapter (Dubiner order). On several processes `local_mat` holds partial
rows (each process's own cells only, no global numbering): the adapter extracts the blocks per process, numbers the
faces (exclusive scan) and sums the rows of the shared faces on their owners (`.scratch/ngsolve/issues/03-...`).
Non-symmetric: one process only. Towards hypre and PETSc NGSolve sends its partial rows and lets them sum
(`HYPRE_IJMatrixAddToValues`; ngsPETSc: MATIS, then AIJ).

**PETSc users.** Each process assembles its own cells; PETSc sends the contributions to rows of other processes. A
block on several processes is an MPIAIJ: two CSRs (the owned and the other columns) and `garray`, the global DoF
numbers of the other columns. A copy before the library is accepted (as for hypre, AmgX, MUMPS). The interpolation of
1 is the near-null space vector (`MatGetNearNullSpace`). fhhos4's preconditioner, `-pc_type fhhos4` (`library/petsc`,
loaded by `-dll_append`), takes a MATNEST of the four blocks; a monolithic AIJ matrix with index sets: not yet.

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
`UTILISATEUR`).

From its sources (gitlab.com/codeaster/src, main of 2026-10-21; checked 2026-11-17):
- **The PETSc matrix** (`bibfor/third_party_interf/apalm*.F90`, `apmam*.F90`): one MPIAIJ (SEQAIJ on one process),
  both triangles (a symmetric MORSE matrix stores half: the conversion writes the transpose), all the stored entries,
  zeros included (`ADD_VALUES`, no `MAT_IGNORE_ZERO_ENTRIES`). PETSc with 32- or 64-bit indices
  (`ASTER_PETSC_64BIT_INDICES`).
- **The entry point**: `SOLVEUR=_F(METHODE="PETSC", PRE_COND="UTILISATEUR", KSP_UTIL=f)`: Code_Aster calls `f(mat)`
  (`bibcxx/Solvers/MatrixToPetsc.cxx`, `create_custom_ksp`; petsc4py builds only, `ASTER_HAVE_PETSC4PY`) with its
  matrix as a petsc4py `Mat`, and solves with the KSP returned. Its test `astest/petsc05b.comm` does what we need: the
  script builds `NUM = NUME_DDL(MODELE=..., CHARGE=...)`, then `f` builds index sets from it (`getPhysicalDOFs()`,
  restricted to the ownership range on several processes) and the blocks (`createSubMatrix()`). Without `KSP_UTIL`,
  `PRE_COND="UTILISATEUR"` only inserts `OPTION_PETSC` and calls `PCSetFromOptions()` (`appcpr.F90`). `OPTION_PETSC`
  also goes to `PetscInitialize()` at the first PETSc solve (`apetsc.F90`), where `-dll_append` loads a library.
- **HHO's unknowns** (`catalo/cataelem`): on the center node of each face (`HHO_FX1..n`, `HHO_FY..`, `HHO_FZ..`)
  and of each cell (`HHO_CX..`, `HHO_CY..`, `HHO_CZ..`); meshes with those nodes (`HEXA27`...). The component names
  give the cell and face index sets (`DOFNumbering.getNodeAndComponentFromDOF()`). Within a node, the components
  follow the physical quantity's order: all of x, then all of y (structure of arrays), where the library wants them
  interleaved (`phi_i e_c` is the DoF `n i + c`): the face index set must list them in the library's order.
  Gradient-enhanced models (`_GRAD_HHO`) add `HHO_VR`, `HHO_LG` unknowns.
- **The bases** are L2-orthonormal by default (`BASIS_ORTHO`, `bibfor/hho/HHO_basis_module.F90`): the interpolation
  of 1 is not ones, so the near null space must be given. `HHO.projectOnHHOSpace()` (Python) projects a function per
  component onto the HHO unknowns: the translations and rotations (`MatSetNearNullSpace`). Code_Aster's own rigid
  body modes for GAMG/ML (`MatNullSpaceCreateRigidBody` from the nodal coordinates) do not fit HHO.
- **Dirichlet conditions**: by elimination (`AFFE_CHAR_CINE`) in 83 of its 88 HHO mechanics tests. The eliminated
  unknowns stay in the matrix: their rows and columns zeroed, 1 on the diagonal (`asmchc.F90`), the zeros stored.
  A Dirichlet face therefore stays a face block, decoupled; a face fixed in one component only is half decoupled.
  `LIAISON_*` couplings and contact add Lagrange multipliers (`ELIM_LAGR="OUI"` exists for PETSc; untried with HHO).
- **Its solvers**: MUMPS by default; with PETSc, `LDLT_SP` (a single-precision MUMPS factorization as preconditioner,
  refactorized when a solve takes more than `REAC_PRECOND` iterations, `pcmump.F90`). `PRE_COND="HPDDM"`: HPDDM with
  GenEO (symmetric; `HARMO` otherwise), only with a distributed mesh (`appcpr.F90`): MUMPS on each subdomain, the coarse
  space from the subdomains' Neumann matrices, assembled from the elementary matrices (`apmams.F90`). The test suite
  uses MUMPS or `LDLT_SP` with HHO and HPDDM only with Lagrange elements, but EDF's studies, not in the repository,
  use HPDDM with HHO (someone at EDF, 2026-11-17): the solver to beat on large cases.

**MFEM, deal.II.** Not studied: from memory.

## What it means for the library

- Blocks assembled by the client (the program, HArDCore3D, Disk++): the setups alone, no helper.
- A monolithic matrix is common (NGSolve, Code_Aster, PETSc's AIJ, DOLFINx's monolithic assembly): `ExtractBlocks*()`
  on one process; on several processes, through PETSc (`MatCreateSubMatrix()`) or, for a client without PETSc, a
  distributed `ExtractBlocks*()` (not planned).
- PETSc's and hypre's distributed matrices (two CSRs per block): `MergeOffDiag()`.
- Partial rows (NGSolve, PETSc's MATIS): summed by the client today.
- Repeated setups on the same pattern (Code_Aster's Newton iterations): the setup time matters.
- Dirichlet unknowns kept in the matrix as decoupled rows (Code_Aster's `AFFE_CHAR_CINE`), possibly for some
  components of a face only: untested in the library.
