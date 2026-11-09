# No new dependencies; framework integrations live outside the repository

fhhos4 takes no new dependency if it can be avoided (2026-10-19: "No new dependency if possible"). An option that needs
one is acceptable only when it is off by default (`ENABLE_PYTHON`: nanobind, NumPy). Tests and benchmarks that need a
client framework (NGSolve, DOLFINx) stay outside the repository, and the repository's documentation does not depend on
them: NGSolve's adapter and its tests live in a directory of their own.

## Bindings to solver libraries (2026-11-09)

The preconditioner for PETSc (`library/petsc`, `-pc_type fhhos4`) lives in the repository, behind `ENABLE_PETSC_PC`
(off): PETSc is a solver library that fhhos4 plugs into, not a framework whose tests need its discretization; the
binding's tests build their system themselves (`tests/PETSc`), its example (`examples/`) its own HHO with petsc4py. A separate library, `libfhhos4_petsc`: `libfhhos4` never
depends on PETSc. Its option is not named `ENABLE_PETSC`, which would read as fhhos4 using PETSc, as it uses MUMPS
(`ENABLE_MUMPS`). Meant to move into PETSc (`src/ksp/pc/impls`, `--download-fhhos4`, like hypre's, AmgX's and HPDDM's
preconditioners) once stable, the file under PETSc's license; then removed from here once the PETSc versions with it
are common.
