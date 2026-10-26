# Mixed precision: the levels in float, their setup in double

With `MixedPrecision` (the library's default, 2026-10-23), the cycle runs in single precision and the outer Krylov
method, residuals and solution in double. The setup is the double one, its results rounded to float as they are
stored, so the aggregates and coarse matrices are those of the double precision. One templated code for both
precisions (`UncondensedLevelOf<T>`): "Do not duplicate the code, use templating".

## Considered Options

- Levels built from float copies of the blocks: the rounded couplings changed the aggregates where couplings nearly
  tie (Square-tri n=16 k=0, the cycles alone: 32 iterations instead of 29). Rejected.
