# Eigen matrices only in the public API, and no C interface

Users see and manage Eigen matrices only (`fhhos4::SparseMatrix` or an Eigen map of one), "even if they lose a little
performance" (2026-10-25); the implementation types stay in `fhhos4::detail`, out of every signature the user calls.
The compiled functions take plain views; the Eigen functions are inline in the header, compiled with the client's
Eigen, so no Eigen object crosses the compiled boundary.

## Considered Options

- A C interface (rejected 2026-10-19): a PETSc preconditioner would be a C++ file, like PETSc's own PCHPDDM and
  PCAMGX; MATNEST blocks are read in place, AIJ plus index sets go through `ExtractBlocksSym()`.
- Raw CSR views in the public API (2026-10-19 to 10-25): kept only for the Python module, behind protected members.
