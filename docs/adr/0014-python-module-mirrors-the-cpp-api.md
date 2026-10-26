# The Python module mirrors the C++ API

The Python module has the C++ names (`solver.FaceDegree`, `SetupSymmetric()`, `Solve()`), so that one documentation
serves both (2026-10-18). It is built in this repository (nanobind, `ENABLE_PYTHON`, off) and reads NumPy arrays in
place like the C++ library (ADR-0003): no conversion, float64 vectors only; SciPy matrices are duck-typed, SciPy is
not needed. Its pip or conda packaging is decided after its validation in the client frameworks.
