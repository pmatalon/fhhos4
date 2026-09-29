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
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

Produces `build/bin/fhhos4` (main solver) and `build/bin/fhhos4_tests` (GTest suite).

## Test

From `build/`:

```bash
ctest --output-on-failure
```

39 tests, ~2 min. Covers `Program_Diffusion_HHO` with `mg`/`fcgmg` solvers on
`cart`/`stri`/`stetra` (in-house) and GMSH `cart` meshes — see `tests/README.md`.
Not covered: biharmonic, DG/FEM programs, other solvers (`lu`, `cg`, `uamg`,
`agmg`, ...), 1D, anisotropic cases.

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
- `Mesher/` — mesh generation: `InHouse/` (structured meshes) and `GMSH/` binding,
  plus built-in geometries (Square, Cube, Square4quadrants).
- `Mesh/` — mesh data structures (Element, Face, Vertex, agglomeration).
- `Discretizations/` — `DG/`, `FEM/`, `HHO/` discretization schemes.
- `Solver/` — `Direct/`, `Krylov/`, `Multigrid/`, `FixedPoint/`, `BiHarmonic/`,
  plus `SolverFactory.h`.
- `FunctionalBasis/`, `QuadratureRules/`, `Geometry/` — numerics building blocks.
- `TestCases/` — analytic/benchmark problem definitions.
- `Utils/` — cross-cutting helpers (parallel loops, timers, export, types).

## Troubleshooting

- `cmake`/`make`/CGAL/Eigen not found → conda env not activated.
- GMP/MPFR not found → `conda install -n fhhos4 gmp mpfr`.
