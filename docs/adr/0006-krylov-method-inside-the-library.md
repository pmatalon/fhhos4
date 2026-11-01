# The Krylov method is inside the library

With the K-cycle, U-AMG is a non-linear preconditioner: paired with CG or another non-flexible method it breaks. So the
library runs its own Krylov method (`Krylov`: FCG by default, GCR, or none; BiCGSTAB, with a V- or W-cycle only, was
removed on 2026-10-30: it was for testing only) and exposes no `Apply()` nor preconditioner adapter (removed
2026-10-06), so that its users cannot fall into that trap.
`Krylov::None` with `MaxIterations = 1` (zero initial guess) is one cycle: that is how a caller preconditions its own
Krylov method, on the condensed system (`CondenseRHS()`, `SolveCellUnknowns()`): a flexible one with the K-cycle; V-
and W-cycles are linear, for any method, but perform poorly.
