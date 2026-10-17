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

Produces `build/bin/fhhos4` (main solver), `build/bin/fhhos4_tests` (GTest suite) and `build/lib/libfhhos4.so`
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
- `ENABLE_PROGRAM` (ON): OFF builds the library `fhhos4` alone (needs only Eigen): the CMake file has a library
  section, then the program's, skipped by `return()`.

`cmake --install . --prefix <prefix> [--component library|program]`: program, meshes (`share/fhhos4/meshes`), library,
headers, CMake package. No path of the sources is compiled in (`Utils/FileSystem.h`): the program finds the meshes
from its own location (`<bin>/../share/fhhos4`, a link to `data/meshes` in the build tree; `FHHOS4_DATA_DIR`
overrides; a data file not found is a `Utils::FatalError(FileSystem::DataFileNotFound(...))`, which tells how to set it), writes its exports in `./out` (working directory, `-o`), its GMSH mesh cache in `~/.cache/fhhos4`
(`XDG_CACHE_HOME`), and its temporary files in the system's temporary directory, named after the process id.

Conda package `fhhos4`: the library alone (`conda/recipe/recipe.yaml`, rattler-build v1 recipe; `rattler-build build -r
conda/recipe -c conda-forge --experimental`, ~1 min; rattler-build is in the conda env `rattler`, not `fhhos4`).
`ENABLE_PROGRAM=OFF`, LGPL. Run dependencies: Eigen >= 3.4 (the user's decision, 2026-10-17: a minimum, no exact pin;
the client's Eigen may differ from the library's, see "U-AMG: details and measurements"), and the compiler's runtime
(`libstdcxx`, `libgcc`; the overlinking warning on `libgomp` is normal on conda-forge). `run_exports` pinned to
major.minor, the SOVERSION. The program is not packaged (the name `fhhos4` is the library's: the clients install it).
Linux only. Its tests: `conda/recipe/test`, a CMake client (version, a solve on a lowest-order hybrid system), compiled
as is, with `-march=native` and with Eigen 3.4. The source is the clone (`path: ../..`, without the gitignored files):
for conda-forge, the archive of the release tag (url and sha256) and the version written out. `variants.yaml`: local
builds only.

## Test

From `build/`:

```bash
ctest --output-on-failure
```

212 tests, ~5 min (`ctest -E HPConfig` skips the 8 slowest). Most check the
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
  plus `SolverFactory.h` (needs the problem), `AlgebraicSolverFactory.h` (only the matrices: shared with the library),
  `StaticCondensation.h` (algebraic, by cells: the library's) and `LibraryUAMG.h/.cpp` (`-s libuamg`).
- `FunctionalBasis/`, `QuadratureRules/`, `Geometry/` — numerics building blocks.
- `TestCases/` — analytic/benchmark problem definitions.
- `Utils/` — cross-cutting helpers (timers, export, types). `Utils::FatalError()` throws `fhhos4::Error` (`main()`
  prints it and exits with EXIT_FAILURE): never call it inside a parallel loop. The solvers print through `Utils::Log()`
  (cout, redirected by the library according to its verbosity). Parallel loops are plain OpenMP
  (`#pragma omp parallel for`); `Parallelism.h` has the thread count (default: all logical cores, but one per
  physical core in the solvers: `SolverThreads`, around the solver section of each program) and `ThreadLocal`/
  `ThreadLocalCoeffs` for per-thread results (e.g. the non-zeros of a matrix being assembled).

## Library (`library/`)

`fhhos4` (shared library `libfhhos4.so`; CMake target `fhhos4_lib`, since `fhhos4` is the program; exported as
`fhhos4::fhhos4`) gives U-AMG to other HHO codes (first target: HArDCore3D), see `README_library.md` (its rules: section
"Documentation"). Public interface in `library/include/fhhos4.h` (`#include <fhhos4>`, a
file that includes `fhhos4.h`, or `<fhhos4.h>`; namespace `fhhos4`, Eigen types, pimpl): `Solver` (parameters as
public members; `Dimension`, `FaceDegree` and `CellDegree` must be set), `Result`, and `Error` (in `fhhos4_error.h`,
which `Utils.h` includes alone). Two interfaces: the uncondensed system, the interface of choice
(`SetupSymmetric()` from `A_TT`, `A_TF`, `A_FF` and the interpolations of 1, or `SetupNonsymmetric()` with `A_FT` too,
then `Solve(b_T, b_F, x_T, x_F)`: the library condenses the system by cells, `StaticCondensation.h`, and computes the
cell unknowns; `CondensedMatrix()`), and the condensed system, only to integrate the library with the least effort into
codes that condense already (HArDCore3D, Disk++: `SetupSymmetricCondensed()`/`SetupNonsymmetricCondensed()` from `A`
and the blocks, `A_FF` optional, then `Solve(b, x)`). `fhhos4::SparseMatrix`: fhhos4's matrix type, its index type
generated into `fhhos4_version.h` (`SMALL_INDEX`). Non-symmetric matrix (`UncondensedLevel::Symmetric` false, `A_F_T` given; program option
`-uamg-nonsym` to test it on the symmetric matrices): the same algorithm on the matrix itself, without the optimizations of the symmetric case (lower triangular
parts, Cholesky; the local Galerkin products compute full matrices); couplings from the averaged traces of `A_TF` and
`A_FT`, Galerkin coarse operators `P^T A P`, Gauss-Seidel sweeps in the downwind order of each level
(`Solver/Multigrid/DownwindOrder.h`, pre-
and post-smoothing), K-cycle with GCR. Implementation `library/Solver.cpp`, compiled with
`-fvisibility=hidden` (only the public API is exported). It builds the solver as the program does
(`AlgebraicSolverFactory`, `ApplyUncondensedAMGDefaults`): its path must not include mesh or discretization headers,
nor read the global `Utils::ProgramArgs`. Its parameters with a choice are enums (`fhhos4::Krylov`, `Cycle`, `Smoother`,
`CoarseSolver`; no strings), as the program's multigrid parameters (`SmootherType`, `CycleType` in `Utils/Types.h`; the
CLI codes in `SmootherFactory`, the coarse solver stays a solver code): `Solver/LibraryParameters.h` maps them both
ways (one table per parameter), for `library/Solver.cpp` and `LibraryUAMG.cpp`. U-AMG is a non-linear preconditioner with the K-cycle: the Krylov method is
inside the library (`Krylov`: FCG by default, GCR, BiCGSTAB with a V/W-cycle only, or none), which gives no cycle
alone. The program runs it with `-s libuamg|fcglibuamg|gcrlibuamg` (`Solver/LibraryUAMG.cpp`, the only file of the
program that includes the public header, through `LibraryParameters.h`: changing it doesn't recompile `Program.cpp`) through the uncondensed
interface (`Diffusion_HHO::CondensationBySolver`: no static condensation in the program), or the condensed one with
`-libuamg-condensed` (and for a full Neumann problem and the exports of the system or of the algebraic error), which
`tests/UAMGLibraryTest.cpp` compares with `uamg|fcguamg|gcruamg`.

### U-AMG: details and measurements (kept out of `README_library.md`)

- Condensation by cells (setup from the blocks): computes the rows of `A` in the pattern of `A_FF`, completed (and
  sorted) only in the rows where `A_FF` lacks couplings of the faces of a cell (zeros dropped by an assembly): 1.3x (1
  thread) to 1.8x (8 threads) slower when all the rows are completed. HArDCore3D 32^3, k = 1 (285,696 face DoFs, 9.2 M
  non-zeros in `A`): setup from the blocks 0.71-0.80 s vs 0.38-0.50 s from `A` on 1 thread, 0.25-0.35 s vs 0.17-0.29 s
  on 8; same solve. The setups from the blocks keep the inverses of the blocks of `A_TT`, `A_TT^{-1} A_TF` and `A_FT`.
- p-levels: `P_CS` Minus2 (k, k-2, ..., 0), blocks extracted (first coefficients of each block,
  `SetupOperatorByBlockExtraction()`), prolongation = injection of the low modes. h-levels: multiple pairwise
  aggregation (factor 3.8), interface collapsing, coarsening prolongation `ReconstructSmoothedTraceOrInject` (trace of
  `-A_TcTc^{-1} A_TcFc` on the removed faces, smoothed by block Jacobi omega = 2/3), chained as multigrid prolongation.
- Interpolation of 1 (one coefficient c per cell and per face, 1 = c phi_0; the same on the p-levels): the AMG rescales
  the constant modes before the h-coarsening (`UncondensedLevel::CoarsenMesh()`): without it, L2-orthonormal bases slow
  it down (Cube-tet n=16, k=0: 44 iterations instead of 26) or make it fail (graded meshes).
- Non-symmetric matrix: couplings from the averages of the traces of `A_TF(T, F)` and `A_FT(F, T)`; the couplings of
  the rows of the matrix (as AGMG's strong couplings, from `A_TF(T1, F)` and `A_FT(F, T2)`) cost 20% more iterations
  on HArDCore3D's diffusion-dominated convection-diffusion at k = 0, more without the downwind order. Reconstruction
  `-A_TT^{-1} A_TF` by LU of the blocks; smoothing of the prolongations with the rows of the matrix; Galerkin `P^T A P`;
  defaults GCR and LU (FCG and Cholesky refused). Downwind order (Bey & Wittum 1997, Hackbusch & Probst 1997): faces
  sorted after those they are strongly and one-sidedly coupled to, cycles broken; natural order without one-sided
  couplings; on a level with one-sided couplings all the Gauss-Seidel sweeps (backward, symmetric included) follow it,
  Jacobi and SOR ignore it. HArDCore3D convection-dominated 32^3, k = 1: 16 iterations, 4 with it.
- GCR: restarted every 10 iterations, Notay's economical version (ETNA 37, 2010, Algorithm 3.1); the K-cycle then
  iterates with GCR. BiCGSTAB is not flexible: the setup refuses the K-cycle with it.
- Default smoothers: block Chebyshev of the fourth kind, degree 2, for a symmetric matrix (below); for a
  non-symmetric one, on several threads, hybrid block Gauss-Seidel (`hbgs`, `hrbgs`: Gauss-Seidel in the rows of each
  thread, Jacobi between the threads; scalar kernels at k = 0), results depending on the number of threads, not on the
  run; sequential sweeps on a level where more than 20% of the non-zeros couple rows of different threads (e.g. GMSH's
  numbering on the tetrahedral meshes). 0-6 more iterations than the block Gauss-Seidel with FCG (a few more with
  BiCGSTAB), but faster on 8 threads: 15-20% at k = 0, 35-45% at k = 1, 2 (HArDCore3D cubes and Voronoi); BiCGSTAB or
  the multigrid alone 2-45% (earlier, non-deterministic version). GCR, HArDCore3D diffusion-dominated
  convection-diffusion 48^3, k = 1: 27 iterations on 4 threads (26 on 1), solve 3.3-3.5 s instead of 5.8 s.
- Chebyshev smoothers (`bcheb4`, `bcheb`; `Solver/FixedPoint/BlockChebyshev.h`): symmetric matrices only (the library
  refuses them otherwise); `lambdaMax` of `D^-1 A` by 10 Lanczos iterations at the setup. Faster than `hbgs` on 8
  threads in the program (up to 2x on GMSH tetrahedra at k >= 1) and in Disk++ (1.1-2x at k >= 1), a tie in
  HArDCore3D (up to 1.2x on Voronoi meshes), slower on 1 thread: `PERFORMANCE.md`. The library's
  default for a symmetric matrix, on any number of threads (`Smoother::BlockChebyshevFourthKind`, `bcheb4` in the
  program; `Pre/PostSmoothingIterations` -1: degree 2; the
  user's decision, 2026-10-15: the same results whatever the number of threads). The program keeps `bgs,rbgs` (the
  papers).
- FCG replaces its recurrence residual by the true `b - A x` every 20 iterations and when it meets the tolerance (GCR:
  at each restart and at the tolerance), then continues if the true one does not: the tolerance applies to the true
  residual. Near the attainable accuracy (~1e-11 relative on Disk++'s 2D systems at k = 2) the recurrence drifts (8e-13
  for a true 1.4e-11); a tolerance below it runs `MaxIterations`, not converged.
- Symmetric condensed interface: FCG applies `A` from its lower triangle; for an `A` symmetric up to rounding (local
  Schur complements), the solution differs from a direct solver's on the full `A` by ~cond(A) x asymmetry (Disk++ 2D,
  k = 2, asymmetry 3e-16: up to 1% of L2 errors of 3e-11).
- `Threads` default: at most one thread per physical core unless `OMP_NUM_THREADS` or OpenMP's binding is set (the
  hyperthreads slow the AMG down). Release builds use the conda compilers' `-march=nocona`.
- Eigen across the interface (2026-10-17, replaces the comparison of the version and of `EIGEN_MAX_ALIGN_BYTES`): the
  constructor compares the sizes of `Eigen::Index`, `SparseMatrix` and the two `Eigen::Ref`, and whether
  `Ref<VectorXd>` (passed by value) is trivially copyable, then reads samples the header builds with the client's Eigen.
  Eigen 3.4.0 and 5.0.1 have the same layouts (members, sizes, mangled names), and their aligned allocators differ,
  harmless since no memory changes hands; nothing reads the client's data as aligned (`Ref` is unaligned by default).
  Eigen 3.4 only lacked `saveMarketDense` (`Utils/ExportModule.h` writes the format itself). Eigen 5's legacy macros
  say 3.5.0 (`EIGEN_WORLD_VERSION` stays 3; `EIGEN_PATCH_VERSION` exists from 5.0).

## Version

Set in `project()` (`CMakeLists.txt`) only; generated from it: the public header `fhhos4_version.h`
(`FHHOS4_VERSION`...), the CMake package version (`find_package(fhhos4 1.4)`: same major version), the library's
SOVERSION (major.minor: adding a public member to `fhhos4::Solver` breaks the binary compatibility, so only in a new
minor or major version) and the conda package's version (the recipe reads `project()` with `load_from_file`, an
experimental feature of rattler-build: `--experimental`). `fhhos4::Version()` (printed by the program at start) adds
`git describe` when built from a clone: `cmake/GitVersion.cmake` regenerates `build/generated/Version.cpp` at each
build, rewritten only when the commit changes; nothing else may depend on the commit (it would recompile `Program.cpp` at each commit). Release: change the
version in `project()`, commit, then push an annotated tag `vX.Y.Z` (the tags before 1.3.0 are `Release_X.Y.Z`):
`.github/workflows/release.yml` checks the tag against the version, builds, runs `ctest`, the installation and the
library alone, and in a parallel job the conda package with its tests, then creates the GitHub release "Release X.Y.Z"
(notes: the message of the tag). It pins GMSH, Eigen, CGAL, the compiler and rattler-build to the versions the tests
were validated with: update them with the local environment (envs `fhhos4` and `rattler`).

## Documentation (README files)

Three READMEs, each for its readers; keep them up to date in the same commit as the change they describe:
- `README.md`: minimal. What fhhos4 is, its two parts (one paragraph each, linking to the two files below), and the
  license. Nothing else.
- `README_program.md`: the program, for its users (researchers running experiments): dependencies, conda environment,
  build (options table: every `option()` of `CMakeLists.txt`), installation, where the program reads and writes,
  troubleshooting.
- `README_library.md`: the library, for the developers of client codes, who don't need the program or the sources:
  the algorithm in a few lines (p-levels, then h-levels and how the blocks give the aggregation), the interfaces,
  install from source (library alone), link (CMake package), use (examples of both interfaces, non-symmetric), inputs,
  parameters, constraints, license, version (last). Rules:
  - no details of the algorithm beyond that, no measurements, no client codes (HArDCore3D, Disk++), no internal
    files (`PERFORMANCE.md`, tests): they go in "U-AMG: details and measurements" above or in `PERFORMANCE.md`;
  - the parameters table gives the library's own defaults (`library/include/fhhos4.h` and `library/Solver.cpp`:
    `Arguments()`, `CheckParameters()`, `DefaultPreSmoother()`), which differ from the program's (smoother on several
    threads, true residual, the defaults of the non-symmetric case): read them in the code, never from the program;
  - any change of the public interface (`library/include/*`: a setup or solve function, a parameter, its default or
    its accepted values, an input requirement, an error) or of the library's build, install or CMake package updates
    `README_library.md` and the comments of `fhhos4.h`, which document the same interface; a version change updates
    its Version section.

## License

LGPL-3.0-or-later (`COPYING.LESSER`, `COPYING`; README section "License"). Copyright: CERFACS for 2018-2022, Pierre
Matalon since 2023. The public headers (`library/include`), installed into other codes' prefixes, start with an
`SPDX-License-Identifier` line and their copyright line; the other files have none. Code taken from elsewhere must be
compatible with the license. The library must not depend on the GPL dependencies of the program (GMSH, CGAL): codes
under other licenses link it.

## Troubleshooting

- `cmake`/`make`/CGAL/Eigen not found → conda env not activated.
- GMP/MPFR not found → `conda install -n fhhos4 gmp mpfr`.
