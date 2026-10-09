# fhhos4

Fast HHO Solver for Research — solves diffusion (and biharmonic) PDEs with Hybrid
High-Order (HHO) and DG discretizations, with a focus on geometric/algebraic
multigrid solvers. Numerical experiments from related papers are reproduced via
`reproducibility/*.md`.

## Environment

Must activate the conda env before building or running anything (provides Eigen,
CGAL, GMSH, GTest, the compiler, cmake, make):

```bash
conda activate fhhos4
```

If the env doesn't exist yet: `conda env create -f conda/environment.yml`.

## Build

```bash
conda activate fhhos4
cd build          # or: mkdir build && cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ..
ninja
```

Produces `build/bin/fhhos4` (main solver), `build/bin/fhhos4_tests` (GTest suite) and `build/lib/libfhhos4_AMG.so`
(U-AMG for other codes, see `library/`).
Clean build ~2.5 min. Ninja runs at most `MAX_COMPILE_JOBS` (default 2) compilations in
parallel: `Program.cpp` alone takes ~4 GB of RAM and ~2.5 min, and sets the build time.
With Make (`cmake ..` without `-G Ninja`), use `make -j2`.

CMake options (`-D<OPTION>=ON/OFF`):
- `ENABLE_2D`, `ENABLE_3D` (ON), `ENABLE_1D` (OFF): dimensions to compile.
- `ENABLE_DG`, `ENABLE_FEM` (OFF): DG and FEM programs (`-discr dg|fem`).
- `ENABLE_TESTS` (ON): build `fhhos4_tests`. Turning it off saves little time with Ninja
  (the tests compile in parallel with `Program.cpp`), more with Make.
- `ENABLE_GMSH`, `ENABLE_CGAL` (ON, CGAL is only used in 2D), `ENABLE_AGMG` (OFF).

## Test

From `build/`:

```bash
ctest --output-on-failure
```

130 tests, ~4 min (`ctest -E HPConfig` skips the 8 slowest). Most check the
iteration counts of the papers in `reproducibility/` — see `tests/README.md`. Covers
`Program_Diffusion_HHO` and `Program_BiHarmonic_HHO` (HHO, static condensation), 2D/3D,
in-house, GMSH and polygonal meshes, and the library (`UAMGLibrary*`, ~15 s). Not covered: DG/FEM programs (disabled by default), `lu`/`cg`/`agmg`/`p_mg` solvers,
1D, anisotropic cases.

## Performance

Performance is key for this project. Before optimizing, read `PERFORMANCE.md`: current profile, how to measure and
validate a change (the papers' iteration counts must stay unchanged; reordering arithmetic operations is allowed if it
helps performance, changing the numerical method must be opt-in), and the ranked list of next ideas with their
expected gains. Update it after measuring or optimizing.

Measured, so don't re-litigate without new data: Eigen's kernels (sparse matrix x vector, vector operations, small
dense block solves) run at the hardware limit, and its row-major sparse matrix x vector product is already parallel
(OpenMP). Hand-written copies or BLAS/LAPACK don't make them faster. The gains come from reading the matrices fewer
times (fused passes), skipping work, storing fewer bytes, and parallelizing what Eigen does sequentially (sparse x
sparse products, transposes).

## Run

```bash
./bin/fhhos4 -h                # full CLI help
./bin/fhhos4 -geo square -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k 0 -n 32
```

Example invocations for specific experiments/figures live in `reproducibility/*.md`
(one file per paper) and `tests/README.md` (one command per test case).

## Source layout (`src/`)

- `Program/` — top-level drivers per problem type (`Program_Diffusion_HHO/DG/FEM`,
  `Program_BiHarmonic_*`); `main.cpp` and `ProgramArguments*.h` are the CLI entry
  point and argument parsing/defaulting.
  `Program/Program_*.h` and `Program.h` only declare the programs: their definitions
  (`Program/Program_*_Impl.h`) are compiled once, in `Program.cpp`, with explicit
  instantiations per enabled dimension. Only include the light headers elsewhere (main,
  tests): including an `_Impl.h` recompiles the whole solver in that translation unit.
- `Mesher/` — mesh generation: `InHouse/` (structured meshes) and `GMSH/` binding,
  plus built-in geometries (Square, Cube, Square4quadrants).
- `Mesh/` — mesh data structures (Element, Face, Vertex, agglomeration).
- `Discretizations/` — `DG/`, `FEM/`, `HHO/` discretization schemes.
- `Solver/` — `Direct/`, `Krylov/`, `Multigrid/`, `FixedPoint/`, `BiHarmonic/`,
  plus `SolverFactory.h` (needs the problem), `AlgebraicSolverFactory.h` (only the matrices: shared with the library)
  and `LibraryUAMG.h/.cpp` (`-s libuamg`).
- `FunctionalBasis/`, `QuadratureRules/`, `Geometry/` — numerics building blocks.
- `TestCases/` — analytic/benchmark problem definitions.
- `Utils/` — cross-cutting helpers (timers, export, types). `Utils::FatalError()` throws `fhhos4::Error` (`main()`
  prints it and exits with EXIT_FAILURE): never call it inside a parallel loop. The solvers print through `Utils::Log()`
  (cout, redirected by the library according to its verbosity). Parallel loops are plain OpenMP
  (`#pragma omp parallel for`); `Parallelism.h` has the thread count (default: all logical cores, but one per
  physical core in the solvers: `SolverThreads`, around the solver section of each program) and `ThreadLocal`/
  `ThreadLocalCoeffs` for per-thread results (e.g. the non-zeros of a matrix being assembled).

## Library (`library/`)

`fhhos4_AMG` (short name fhhos4, the name client codes use): shared library giving U-AMG to other HHO codes (first
target: HArDCore3D), see `library/README.md`. Public interface in `library/include/fhhos4` (`#include <fhhos4>`,
namespace `fhhos4`, Eigen types, pimpl): `Solver` (parameters as public members; `Dimension` and `FaceDegree` must be
set; `Setup()` from `A`, `A_TT`, `A_TF` (`A_FF` optional) and the interpolations of 1, `Solve()`), `Result`, and `Error`
(in `fhhos4_error.h`, which `Utils.h` includes alone). Implementation `library/Solver.cpp`, compiled with
`-fvisibility=hidden` (only the public API is exported). It builds the solver as the program does
(`AlgebraicSolverFactory`, `ApplyUncondensedAMGDefaults`): its path must not include mesh or discretization headers,
nor read the global `Utils::ProgramArgs`. U-AMG is a non-linear preconditioner with the K-cycle: the Krylov method is
inside the library (`Krylov`: FCG by default, BiCGSTAB with a V/W-cycle only, or none), which gives no cycle alone.
The program runs it with `-s libuamg|fcglibuamg` (`Solver/LibraryUAMG.cpp`, the only file of the program that includes
the public header: changing it doesn't recompile `Program.cpp`), which `tests/UAMGLibraryTest.cpp` compares with
`uamg|fcguamg`.

## License

LGPL-3.0-or-later (`COPYING.LESSER`, `COPYING`; README section "License"). Copyright: CERFACS for 2018-2022, Pierre
Matalon since 2023. The public headers (`library/include`), installed into other codes' prefixes, start with an
`SPDX-License-Identifier` line and their copyright line; the other files have none. Code taken from elsewhere must be
compatible with the license. The library must not depend on the GPL dependencies of the program (GMSH, CGAL): codes
under other licenses link it.

## Troubleshooting

- `cmake`/`make`/CGAL/Eigen not found → conda env not activated.
- GMP/MPFR not found → `conda install -n fhhos4 gmp mpfr`.
