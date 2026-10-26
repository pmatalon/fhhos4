# Static condensation by factorizations, never inverses

The library's static condensation keeps the factorizations of the blocks of `A_TT` (Cholesky, or LU with partial
pivoting for a non-symmetric system), not their inverses (2026-10-23): the rounding of explicit inverses limits the
attainable accuracy of the discretization (1e-11..1e-12 instead of 1e-13..1e-14 against a direct solver at k = 3-4).
The smoothers may apply inverses of their diagonal blocks: their rounding does not reach the solution.
