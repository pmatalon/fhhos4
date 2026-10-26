# The Krylov method is inside the library

With the K-cycle, U-AMG is a non-linear preconditioner: paired with CG or another non-flexible method it breaks. So the
library runs its own Krylov method (`Krylov`: FCG by default, GCR, BiCGSTAB with a V- or W-cycle only, or none) and
exposes no `Apply()` nor preconditioner adapter (removed 2026-10-06), so that its users cannot fall into that trap.
`Krylov::None` with one V- or W-cycle is linear: that is how a caller composes it with its own method.
