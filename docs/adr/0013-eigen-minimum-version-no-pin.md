# Eigen >= 3.4 as a minimum, never an exact pin

The packages require Eigen >= 3.4 and accept any later version, also one that differs from the library's (2026-10-17).
This holds because no Eigen object crosses the compiled boundary (ADR-0004; it replaced a check of the layout of the
Eigen objects at run time, 2026-10-19). One dependency on Eigen's internals remains: `BorrowedSparseMatrix` points the
protected members of Eigen's `SparseMatrix` and `CompressedStorage` to the caller's arrays; they are the same in 3.4.0
and 5.0.1, which `tests/BorrowedSparseMatrixTest.cpp` checks. Check it again with each new major Eigen.
