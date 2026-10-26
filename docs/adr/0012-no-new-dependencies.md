# No new dependencies; framework integrations live outside the repository

fhhos4 takes no new dependency if it can be avoided (2026-10-19: "No new dependency if possible"). An option that needs
one is acceptable only when it is off by default (`ENABLE_PYTHON`: nanobind, NumPy). Tests and benchmarks that need a
client framework (NGSolve, DOLFINx) stay outside the repository, and the repository's documentation does not depend on
them: NGSolve's adapter and its tests live in a directory of their own.
