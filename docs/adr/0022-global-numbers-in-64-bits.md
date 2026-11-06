# Global numbers in 64 bits whatever SMALL_INDEX

The global numbers of a distributed system are `GlobalNumber`, `int64_t`, whatever `SMALL_INDEX` (2026-11-04, the
user). This covers the offsets and ghosts of `IndexMap`, the global columns of the rows exchanged between processes or
given to MUMPS (`GlobalColumnsMatrix`), and the messages that carry them. The local numbers (`BigNumber`,
`SparseMatrixIndex`) stay `int` under `SMALL_INDEX`, the default: each process numbers fewer than 2^31 entities, but all
the processes together may number more. Before, the global numbers were `BigNumber`. Beyond 2^31 - 1 faces in total,
the offsets overflowed silently, and the library narrowed its callers' `int64` ghost numbers to `int`, so a ghost of
2^32 + 5 was accepted as the face 5.

## Considered Options

- `BigNumber` renamed `LocalBigNumber`, with `GlobalBigNumber` beside it (the user's first proposal): `BigNumber` has
  some 1300 uses, nearly all local. `GlobalNumber` also matches the public API's term (`faceGhostGlobalNumbers`).
- The global numbers in the library's index type (`SparseIndex`): the client's ghost arrays would change type with the
  library's build, and the default build would cap a distributed system at 2^31 faces, the very case MPI is for.

## Consequences

- The solve loops (products, smoothers, condensation) keep their local 32-bit indices. The matrices with global columns
  exist only during the setup, where their indices take twice the memory.
- MUMPS's indices are `MUMPS_INT`, 32 bits unless MUMPS is built with 64-bit integers: the coarsest level must have
  fewer than 2^31 rows (checked).
- A process still holds fewer than 2^31 DoFs and non-zeros under `SMALL_INDEX`.
- The GMSH tags of the program's meshes (vertices, elements) are exchanged as `GlobalNumber`, but still looked up as
  `BigNumber`.
