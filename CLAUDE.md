# fhhos4

Fast HHO Solver for Research — solves diffusion (and biharmonic) PDEs with Hybrid
High-Order (HHO) and DG discretizations, with a focus on geometric/algebraic
multigrid solvers. Numerical experiments from related papers are reproduced via
`reproducibility/*.md`. Two parts: the program (`src/`, the research test bed) and the library (`library/`, U-AMG for
other codes).

## Where to read before working

Terms: `GLOSSARY.md` (use its words). Decisions: `docs/adr/` (read those of the area before designing a change; flag a
contradiction instead of overriding one). Open work: `.scratch/` (see "Agent skills"). Then, by area:
- **U-AMG**, the static condensation, the smoothers, Krylov methods: `docs/uamg.md`.
- **The library's code** (public header, `library/Solver.cpp`, `-s libuamg`): `docs/library.md`.
- **MPI** (anything distributed): `docs/mpi.md`.
- **The Python module**, and what client frameworks need from it: `docs/python.md`.
- **Integrating the library into a client code**: `docs/library.md`.
- **A release**, the conda recipe, conda-forge, Spack: `docs/packaging.md`.
- **Optimizing or measuring**: `PERFORMANCE.md` first.

## Environment

Must activate the conda env before building or running anything (provides Eigen,
CGAL, GMSH, GTest, the compiler, cmake, make):

```bash
conda activate fhhos4
```

If the env doesn't exist yet: `conda env create -f conda/environment.yml` (MPI and MUMPS included). Build directories:
MPI `docs/mpi.md`; Python, its env and build directories: `docs/python.md`.

## Build

```bash
conda activate fhhos4
cd build          # or: mkdir build && cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ..
ninja
```

Produces `build/bin/fhhos4` (main solver), `build/bin/fhhos4_tests` (GTest suite) and `build/lib/libfhhos4.so`
(U-AMG for other codes, see `library/`).
Clean build ~2.5 min. Ninja runs at most `MAX_COMPILE_JOBS` (default 2) compilations in
parallel: `Program.cpp` alone takes ~4 GB of RAM and ~2.5 min, and sets the build time (ADR-0015).
With Make (`cmake ..` without `-G Ninja`), use `make -j2`.

CMake options (`-D<OPTION>=ON/OFF`):
- `ENABLE_2D`, `ENABLE_3D` (ON), `ENABLE_1D` (OFF): dimensions to compile.
- `ENABLE_DG`, `ENABLE_FEM` (OFF): DG and FEM programs (`-discr dg|fem`).
- `ENABLE_TESTS` (ON): build `fhhos4_tests`. Turning it off saves little time with Ninja
  (the tests compile in parallel with `Program.cpp`), more with Make.
- `ENABLE_GMSH`, `ENABLE_CGAL` (ON, CGAL is only used in 2D), `ENABLE_AGMG` (OFF).
- `ENABLE_PROGRAM` (ON): OFF builds the library `fhhos4` alone (needs only Eigen): the CMake file has a library
  section, then the program's, skipped by `return()`.
- `ENABLE_PYTHON` (OFF): the Python module of the library (`library/python`).
- `ENABLE_MPI` (OFF): the library and the program on several processes; needs MUMPS (the coarse solver, ADR-0021).
- `ENABLE_MUMPS` (OFF; always with `ENABLE_MPI`): MUMPS on one process too (`CoarseSolver::MUMPS`, `mumpsch|mumpslu`).

`cmake --install . --prefix <prefix> [--component library|program]`: program, meshes (`share/fhhos4/meshes`), library,
headers, CMake package. No path of the sources is compiled in (`Utils/FileSystem.h`): the program finds the meshes
from its own location (`<bin>/../share/fhhos4`, a link to `data/meshes` in the build tree; `FHHOS4_DATA_DIR`
overrides; a data file not found is a `Utils::FatalError(FileSystem::DataFileNotFound(...))`, which tells how to set
it), writes its exports in `./out` (working directory, `-o`), its GMSH mesh cache in `~/.cache/fhhos4`
(`XDG_CACHE_HOME`, `-cache-dir`; on several processes also the partitions of the meshes), and its temporary files in the system's temporary directory, named after the process id.

## Test

From `build/`:

```bash
ctest --output-on-failure
```

335 tests, ~8 min (`ctest -E HPConfig` skips the 8 slowest). Most check the
iteration counts of the papers in `reproducibility/` — see `tests/README.md`. Covers
`Program_Diffusion_HHO` and `Program_BiHarmonic_HHO` (HHO, static condensation), 2D/3D,
in-house, GMSH and polygonal meshes, the full Neumann problem (`FullNeumann*`, 2D), and the library (`UAMGLibrary*`,
~15 s). Not covered: DG/FEM programs (disabled by default), `lu`/`cg`/`agmg`/`p_mg` solvers,
1D, anisotropic cases, `SMALL_INDEX=OFF` (64-bit indices; the release workflow only builds the library with it and
runs the conda test client: after a change of index or state handling, build it apart with `-Wtype-limits` and run
the tests there). MPI tests: `docs/mpi.md`; Python: `docs/python.md`.

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
  plus `SolverFactory.h` (needs the problem), `AlgebraicSolverFactory.h` (only the matrices: shared with the library),
  `StaticCondensation.h` (algebraic, by cells: the library's) and `LibraryUAMG.h/.cpp` (`-s libuamg`).
- `FunctionalBasis/`, `QuadratureRules/`, `Geometry/` — numerics building blocks. The tensor bases (monomials,
  Legendre) span Q^q + P^p (`TensorPolynomial.h`): P^k, Q^k (`-poly-space q`) ordered by largest partial degree (the
  first (d+1)^dim span Q^d: hierarchical in Q, what the p-levels of U-AMG need; P^k keeps its order), and their sum.
  HHO in Q (ADR-0016): the reconstruction in Q^l + P^{k+1} (`FunctionalBasisFactory::CreateQPlusP()`: contains the cell
  space, normal derivatives in Q^k(F), Laplacians in Q^l), diffusion only, Cartesian meshes, `-kc 0|-1`, monomials or
  Legendre, no geometric p-multigrid; the orders of P, exact on Q^2 solutions at k = 2 (`HHOConvergenceTest`).
- `TestCases/` — analytic/benchmark problem definitions. The 3D default (sine) test case on `cart -n 4` has zero
  traces on all faces (the condensed right-hand side is rounding noise): compare solvers on n=6 or 8.
- `Utils/` — cross-cutting helpers (timers, export, types). `Utils::FatalError()` throws `fhhos4::Error` (`main()`
  prints it and exits with EXIT_FAILURE): never call it inside a parallel loop. The solvers print through `Utils::Log()`
  (cout, redirected by the library according to its verbosity). Parallel loops are plain OpenMP
  (`#pragma omp parallel for`); `Parallelism.h` has the thread count (default: all logical cores, but one per
  physical core in the solvers: `SolverThreads`, around the solver section of each program) and `ThreadLocal`/
  `ThreadLocalCoeffs` for per-thread results (e.g. the non-zeros of a matrix being assembled).

## Library (`library/`)

`libfhhos4` gives U-AMG to other HHO codes (`README_library.md` for its users; its code: `docs/library.md`). Public
interface `library/include/fhhos4.h`: `Solver` (parameters as public members), `SetupSymmetric()`/
`SetupNonsymmetric()` from the blocks of the uncondensed system (ADR-0002), then `Solve(b_T, b_F, x_T, x_F)`. Inputs
are Eigen matrices only (ADR-0004), read in place, never copied (ADR-0003); inputs whose wrong value silently breaks
the solver have no default (ADR-0005); the Krylov method is inside (ADR-0006); the library's defaults are not the
program's (ADR-0008). Its path must not include mesh or discretization headers, nor read `Utils::ProgramArgs`
(ADR-0001).

## Version

Set in `project()` (`CMakeLists.txt`) only; generated from it: the public header `fhhos4_version.h`
(`FHHOS4_VERSION`...), the CMake package version (`find_package(fhhos4 1.5)`: same major version), the library's
SOVERSION (major.minor: adding a public member to `fhhos4::Solver` breaks the binary compatibility, so only in a new
minor or major version) and the conda package's version. `fhhos4::Version()` (printed by the program at start) adds
`git describe` when built from a clone: `cmake/GitVersion.cmake` regenerates `build/generated/Version.cpp` at each
build, rewritten only when the commit changes; nothing else may depend on the commit (it would recompile
`Program.cpp` at each commit). Releasing: `docs/packaging.md`.

## Documentation

Keep each document up to date in the same commit as the change it describes. Where things go:
- `README.md`: minimal. What fhhos4 is, its two parts (one paragraph each, linking to the two files below), and the
  license. Nothing else.
- `README_program.md`: the program, for its users (researchers running experiments): dependencies, conda environment,
  build (options table: every `option()` of `CMakeLists.txt`), installation, where the program reads and writes,
  troubleshooting.
- `README_library.md`: the library, for the developers of client codes, who don't need the program or the sources. Its
  numbered sections (the user's structure): 1. Installation (install from source, the library alone; the conda
  package; link, CMake package), 2. The algorithm in a nutshell (in a few lines: p-levels, then h-levels and how the
  blocks give the aggregation; applicability), 3. How to use (the uncondensed system, with an example; non-symmetric;
  several processes; Python; inputs, parameters, constraints, license, version last). Rules:
  - no details of the algorithm beyond that, no measurements, no client codes (HArDCore3D, Disk++), no internal
    files (`PERFORMANCE.md`, tests, `docs/`): they go in `docs/*.md` or `PERFORMANCE.md`;
  - the parameters table gives the library's own defaults (`library/include/fhhos4.h` and `library/Solver.cpp`:
    `Arguments()`, `CheckParameters()`, `DefaultPreSmoother()`), read in the code, never from the program (ADR-0008);
  - any change of the public interface (`library/include/*`: a setup or solve function, a parameter, its default or
    its accepted values, an input requirement, an error) or of the library's build, install or CMake package updates
    `README_library.md` and the comments of `fhhos4.h`, which document the same interface; a version change updates
    its Version section.
- `docs/*.md`: how the code works today, for developers and agents; no history, no measurements.
- `GLOSSARY.md`: terms only. `docs/adr/`: a decision that is hard to reverse, surprising, and a real trade-off.
- `PERFORMANCE.md`: every measurement, with its date.
- `.scratch/`: open work (local, gitignored), never in the documents above.

## License

LGPL-3.0-or-later (`COPYING.LESSER`, `COPYING`; README section "License"; ADR-0001). Copyright: CERFACS for 2018-2022,
Pierre Matalon since 2023. The public headers (`library/include`), installed into other codes' prefixes, start with an
`SPDX-License-Identifier` line and their copyright line; the other files have none. Code taken from elsewhere must be
compatible with the license. The library must not depend on the GPL dependencies of the program (GMSH, CGAL).

## Troubleshooting

- `cmake`/`make`/CGAL/Eigen not found → conda env not activated.
- GMP/MPFR not found → `conda install -n fhhos4 gmp mpfr`.

## Agent skills

### Issue tracker

Local markdown: issues and specs are files under `.scratch/<feature>/` (gitignored). See `docs/agents/issue-tracker.md`.

### Triage labels

Default vocabulary (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`), recorded as a
`Status:` line in each issue file. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: one `GLOSSARY.md` and `docs/adr/` at the repo root. See `docs/agents/domain.md`.
