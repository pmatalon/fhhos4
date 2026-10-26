# The solver is compiled in one translation unit

The definitions of the programs (`Program/Program_*_Impl.h`) are compiled once, in `src/Program.cpp`, with explicit
instantiations per enabled dimension (2026-10-01); `main.cpp` and the tests include only the light declarations. Before,
every test file compiled the whole solver: a full rebuild took ~10 min and ran out of memory with parallel jobs. The
price is one heavy file (~4 GB of RAM, ~2.5 min) that sets the build time.

## Consequences

Including an `_Impl.h` anywhere else recompiles the solver there. The program reaches the library's public header
through `Solver/LibraryUAMG.cpp` only, so that changing the header does not recompile `Program.cpp`.
