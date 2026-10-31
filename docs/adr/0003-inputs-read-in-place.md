# The library reads its inputs in place and never copies them silently

The setups and solves read the caller's matrices and vectors in place, in one native format, and never convert or copy
them (2026-10-19). A wrong format is an error that names the conversion (`fhhos4::SparseMatrix B = A;`). Every
conversion is a separate helper documented as a copy (`ExtractBlocksSym()`, `Gather()`, `Scatter()`), so the
caller sees the cost and can change its assembly to avoid it: "keep all copies out of the library". The matrices the
library computes (condensed matrix, levels) are not copies.

## Consequences

The caller's matrices must outlive the solver where it keeps them (`A_FT` for a non-symmetric setup): deleted
overloads refuse temporaries there.
