# One library interface: the uncondensed system

The library takes the blocks of the uncondensed system (`SetupSymmetric()`, `SetupNonsymmetric()`) and condenses it
itself, cell by cell; `Solve(b_T, b_F, x_T, x_F)` returns the cell unknowns too. U-AMG needs the blocks anyway (the
couplings and the coarse levels come from `A_TF`), so the condensed matrix alone is not enough.

## Considered Options

- Setups from the caller's condensed matrix plus the blocks (`Setup...Condensed()`, for codes that condense already):
  kept from 2026-10-14 for the least-effort integrations, removed 2026-10-25 once the clients moved to the uncondensed
  system: one interface to maintain. `Solve(b, x)` on the library's own `CondensedMatrix()` stays.
- Taking `A_TT^{-1} A_TF` instead of `A_TF` (rejected 2026-10-07: "reconstructing A_TF makes little sense").
- Local matrices and a cell-to-face map (rejected 2026-10-07: it needs global indices, not algebraic anymore).
