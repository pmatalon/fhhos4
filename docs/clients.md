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
| PETSc users | the preconditioner (`library/petsc`) | PETSc (`MatSetValues`, `MatAssemblyEnd`) | MATNEST of AIJ blocks, or monolithic AIJ + index sets | yes (or one process) | `-pc_type fhhos4` (`ENABLE_PETSC_PC`: `MergeOffDiag()` from `MatMPIAIJGetSeqAIJ()`); for a monolithic matrix, `PCFHHOS4SetIS()`: `MatCreateSubMatrix()` then the merge |
| DOLFINx | target (`mpi/04`) | itself (`MatrixCSR`) or PETSc (AIJ, MATNEST) | blocks or monolithic matrix | always | through PETSc: as PETSc users; native `MatrixCSR`: to study |
| Code_Aster | target: a prototype (an adapter outside the repository), one process or a distributed mesh | itself (`MATR_ASSE`, MORSE storage); PETSc in its parallel versions | monolithic matrix, uncondensed | yes | through PETSc: a `KSP_UTIL` function of the user's script (petsc4py) builds the index sets and `-pc_type fhhos4` (see below) |
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
loaded by `-dll_append`), takes a MATNEST of the four blocks, or a monolithic AIJ matrix with the index sets of the
cells and faces (`PCFHHOS4SetIS()`: the blocks extracted by `MatCreateSubMatrix()`, the rows left out the identity).

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

The adapter (2026-11-17; conda-forge's `code-aster` 18.1.8, OpenMPI, PETSc 3.25 and petsc4py: build directory
`build-aster`, env `code-aster`): `KSP_UTIL` returns a `preonly` KSP with `-pc_type fhhos4` on Code_Aster's matrix and
the index sets of the cells and faces of each process in the library's order (`PCFHHOS4SetIS()`, composed by
petsc4py), with the library's own FCG to `RESI_RELA` (on the condensed residual); `-ksp_type fcg -pc_fhhos4_krylov
none` in `OPTION_PETSC` for PETSc's FCG around one cycle instead. One process, or several with a distributed mesh
(`ParallelMesh`: `readMedFile()`, `buildCube()`, then `refine()`). The library registered by
`PetscDLLibraryRegister_fhhos4_petsc()` (ctypes). What it relies on:
- The DoF names give the order: `getNodeAndComponentFromDOF()` without argument (all the DoFs in one call; one call per
  DoF takes 260 us), the components interleaved by sorting each node's DoFs by (basis function, component). On a
  distributed mesh (`ParallelDOFNumbering`) these are the local DoFs, ghosts included: `getNoGhostDOFs()` gives the
  owned ones, all those of the nodes the process owns (whole cells and faces), and `getLocalToGlobalMapping()` their
  PETSc rows (0-based, contiguous on each process: `NUME.NULG`, `apalmh.F90`). A sequential mesh on several processes:
  refused by the adapter (untried). The rigid body modes from `HHO(model).projectOnHHOSpace([FORMULE...])` are in
  the matrix's numbering (A v = 0 to rounding away from the Dirichlet rows; on a distributed mesh, on the local DoFs,
  the rotations about a centre reduced over the processes); the bases pass the library's checks (hierarchical,
  constant first), and the discretization is inferred (`-pc_fhhos4_vector_components` given: 3).
- Code_Aster does not set `MAT_SYMMETRIC`, and on several processes its matrix is symmetric to rounding only (the
  beam refined once on 2 processes: |A - A^T|_inf / |A|_inf = 9e-17, on entries that are zero but for rounding): the
  adapter tests it to 1e-12 |A|_inf.
- The eliminated faces (`AFFE_CHAR_CINE`: rows and columns zeroed, 1 on the diagonal) are refused by the library ("a
  face belongs to no cell"): the adapter leaves them out of the blocks, as the library's contract wants the Dirichlet
  faces, and the preconditioner is the identity there. A face fixed in some components only: refused (faces of
  different sizes, `uamg/04`). The stored zeros change nothing (the same residuals without them).
- Code_Aster calls `KSP_UTIL` at each new matrix, keeps the KSP's type, inserts `OPTION_PETSC` and calls
  `PCSetFromOptions()` then `PCSetUp()` (`appcpr.F90`), and only then sets the tolerances (`RESI_RELA`, `NMAX_ITER`:
  10000 if 0) and calls `KSPSetFromOptions()` (`apksp.F90`): the adapter sets the library's Krylov method in
  `KSP_UTIL` (the setup builds it) and passes the tolerances in a pre-solve hook (`KSPSetPreSolve()`). Not the norm type: with `-ksp_type fcg`, `-ksp_norm_type unpreconditioned`, else
  `RESI_RELA` applies to FCG's preconditioned residual (the beam `ssnv303e`: 1e-6 of it is 2e-2 of the true one). After
  the solve, Code_Aster checks the true residual against `sqrt(RESI_RELA)`, in its linear operators only (hidden
  keyword `RESI_RELA_PC`, -1 in `STAT_NON_LINE`: `c_solveur.py`; `apmain.F90`, error `PETSC_16`). With `preonly`, a
  solve of the library that does not converge is not reported (PETSc's `CONVERGED_ITS`).
- Code_Aster sets OpenBLAS's threads itself (`openblas_set_num_threads()` in `libbibfor`, to `ncpus`): on several
  processes the library's MUMPS (one OpenMP thread) would hang in OpenBLAS's products; the library sets OpenBLAS to
  one thread around MUMPS's calls (`docs/mpi.md`).
- A Python object that the KSP keeps (a `PCPYTHON` context, a pre-solve hook) must not refer to it: KSP -> object ->
  KSP is a cycle through C that Python's collector cannot see, and each matrix's KSP and hierarchy would never be
  freed.
- Code_Aster's `PRE_COND="GAMG"` fails on HHO: its rigid body modes (`MatNullSpaceCreateRigidBody()`) need a constant
  number of DoFs per node (`PETSC_18`, a warning; then 10000 iterations). `projectOnHHOSpace()` of a nodal field does
  not exist for mechanics in 18.1.8 (option `HHO_PROJ2_MECA`): the modes from functions, 6 projections per model
  (`PERFORMANCE.md`). conda-forge's PETSc has GAMG and hypre, not HPDDM.
- `create_custom_ksp()` does not check the result of the function: a Python exception there is a segmentation fault
  without its message (the adapter must print it). Code_Aster operators run inside `KSP_UTIL` (`projectOnHHOSpace()`)
  invalidate the matrix it was given ("Not for factored matrix"): compute before the solve.
- PETSc is initialized by the import of `code_aster.Utilities.PETSc` (petsc4py), before `OPTION_PETSC` can reach
  `PetscInitialize()`: the initialization options through the environment (`PETSC_OPTIONS=-no_signal_handler`).

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
