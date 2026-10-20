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
- `ENABLE_PYTHON` (OFF): the Python module of the library (`library/python`, see "Python module" below).

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
as is, with `-march=native` and with Eigen 3.4. The source is the clone (`path: ../..`, without the gitignored files).
conda-forge: staged-recipes PR #35238 (2026-10-17), then the feedstock `conda-forge/fhhos4-feedstock` holds its own
recipe, this one with the version written out and the archive of the release tag as source (url and sha256; its bot
updates them at each release). Both take the test client from the sources (`files: source: conda/recipe/test/`): its
changes reach conda-forge with the release; any other change of `recipe.yaml` goes to the feedstock too.
`variants.yaml`: local builds only.

Spack package `fhhos4`: the library alone too, in Spack's builtin repository (`spack/spack-packages`,
`repos/spack_repo/builtin/packages/fhhos4/package.py`), not in this repository. `CMakePackage`, `ENABLE_PROGRAM=OFF`,
`ENABLE_NATIVE_ARCH=OFF` (Spack sets the target's flags), variant `openmp` (default on), `eigen@3.4:`. No variant for
the index type: `SMALL_INDEX=OFF` does not compile in 1.4.0 (fixed after it): add one (`int64`, `when=` the
next version) with the next version. Its stand-alone test (`spack test run fhhos4`) builds and
runs the conda test client `conda/recipe/test`, cached from the sources: moving it changes the package too (per version).

## Test

From `build/`:

```bash
ctest --output-on-failure
```

221 tests, ~5 min (`ctest -E HPConfig` skips the 8 slowest). Most check the
iteration counts of the papers in `reproducibility/` — see `tests/README.md`. Covers
`Program_Diffusion_HHO` and `Program_BiHarmonic_HHO` (HHO, static condensation), 2D/3D,
in-house, GMSH and polygonal meshes, and the library (`UAMGLibrary*`, ~15 s). Not covered: DG/FEM programs (disabled by default), `lu`/`cg`/`agmg`/`p_mg` solvers,
1D, anisotropic cases, `SMALL_INDEX=OFF` (64-bit indices; the release workflow only builds the library with it and
runs the conda test client: after a change of index or state handling, build it apart with `-Wtype-limits` and run
the tests there).

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
file that includes `fhhos4.h`, or `<fhhos4.h>`; namespace `fhhos4`, pimpl): `Solver` (parameters as
public members; `Dimension`, `FaceDegree` and `CellDegree` must be set), `Result`, and `Error` (in `fhhos4_error.h`,
which `Utils.h` includes alone). Two interfaces: the uncondensed system, the interface of choice
(`SetupSymmetric()` from `A_TT`, `A_TF`, `A_FF` and the interpolations of 1, or `SetupNonsymmetric()` with `A_FT` too,
then `Solve(b_T, b_F, x_T, x_F)`: the library condenses the system by cells, `StaticCondensation.h`, and computes the
cell unknowns; `CondensedMatrix()`), and the condensed system, only to integrate the library with the least effort into
codes that condense already (HArDCore3D, Disk++: `SetupSymmetricCondensed()`/`SetupNonsymmetricCondensed()` from `A`
and the blocks, `A_FF` optional, then `Solve(b, x)`). Inputs (the user's decisions, 2026-10-19): read in place, never
copied nor converted by the setups and solves: the matrices as `CSRView` (CSR arrays of the library's index type
`SparseIndex`, generated into `fhhos4_version.h` from `SMALL_INDEX`; `fhhos4::SparseMatrix`, the Eigen row-major type,
converts to it; checked by one parallel pass, `CheckCSR()`), which the solvers read as `BorrowedSparseMatrix`
(`Utils/BorrowedSparseMatrix.h`: a `SparseMatrix` pointing Eigen's protected members to the caller's arrays); the
library keeps `A` (condensed setups) and `A_FT` (`SetupNonsymmetric()`: `StaticCondensation` reads it at each solve).
Other formats through helpers that copy, so that the copies stay visible in the calling code: `ExtractSymmetricBlocks()`/
`ExtractNonsymmetricBlocks()` (a global matrix + `DoFLists` of the cells and faces, any index types, one parallel pass:
`SparseMatrixOps::BuildByRows`), `Gather()`/`Scatter()`, `ToCSR()` (index types; column-major: the parallel
`SparseMatrixOps::Transpose`); `TransposedView()` (no copy: a symmetric column-major matrix as its transpose). Their
results: `CSRMatrix`, owned by the library. The compiled functions take plain types only (views, `detail::` vector
views): the Eigen overloads and conversions are inline in the header, compiled with the client's Eigen; no C interface
(the user's decision: a PETSc preconditioner would be C++, as PETSc's PCHPDDM and PCAMGX; MATNEST blocks read in place,
AIJ + index sets through `ExtractSymmetricBlocks()`). Non-symmetric matrix (`UncondensedLevel::Symmetric` false, `A_F_T` given; program option
`-uamg-nonsym` to test it on the symmetric matrices): the same algorithm on the matrix itself, without the optimizations of the symmetric case (lower triangular
parts, Cholesky; the local Galerkin products compute full matrices); couplings from the averaged traces of `A_TF` and
`A_FT`, Galerkin coarse operators `P^T A P`, Gauss-Seidel sweeps in the downwind order of each level
(`Solver/Multigrid/DownwindOrder.h`, pre-
and post-smoothing), K-cycle with GCR. Implementation `library/Solver.cpp`, compiled with
`-fvisibility=hidden` (only the public API is exported). It builds the solver as the program does
(`AlgebraicSolverFactory`, `ApplyUncondensedAMGDefaults`): its path must not include mesh or discretization headers,
nor read the global `Utils::ProgramArgs`. Its parameters with a choice are enums (`fhhos4::Krylov`, `Cycle`, `PCoarsening`,
`Smoother`, `CoarseSolver`; no strings), as the program's multigrid parameters (`SmootherType`, `CycleType`, `P_CoarsStgy` in `Utils/Types.h`; the
CLI codes in `SmootherFactory`, the coarse solver stays a solver code): `Solver/LibraryParameters.h` maps them both
ways (one table per parameter), for `library/Solver.cpp` and `LibraryUAMG.cpp`. U-AMG is a non-linear preconditioner with the K-cycle: the Krylov method is
inside the library (`Krylov`: FCG by default, GCR, BiCGSTAB with a V/W-cycle only, or none), which gives no cycle
alone. The program runs it with `-s libuamg|fcglibuamg|gcrlibuamg` (`Solver/LibraryUAMG.cpp`, the only file of the
program that includes the public header, through `LibraryParameters.h`: changing it doesn't recompile `Program.cpp`) through the uncondensed
interface (`Diffusion_HHO::CondensationBySolver`: no static condensation in the program), or the condensed one with
`-libuamg-condensed` (and for a full Neumann problem and the exports of the system or of the algebraic error), which
`tests/UAMGLibraryTest.cpp` compares with `uamg|fcguamg|gcruamg`.

### Python module

`library/python/fhhos4_python.cpp`: the module `fhhos4` (nanobind, `NB_STATIC`), the public API with the C++ names (the
user's decision, 2026-10-18: one documentation), on NumPy arrays read in place, as the library reads its inputs (never
copied nor converted: `nb::try_cast` without conversion; the vectors `noconvert()`, float64 only): the matrices as
`CSRView` (`PyCSRView`: the arrays of any index types, int32/int64, uint64 read as int64, and their Python owners),
`CSRMatrix` (arrays exposed as read-only NumPy views that keep it alive), or any object with `indptr`, `indices`,
`data`, `shape` and `format` "csr" (SciPy, duck-typed: SciPy is not needed, only by the test). `PySolver` derives from
`fhhos4::Solver` and keeps the views the solves read (`A`, `A_FT`), replaced only after a successful setup: the checks
of a setup run before it resets the previous one, which still reads them. `CondensedMatrix()`: a copy (`ToCSR()`). The
setups, solves and helpers release the GIL. Built with `-DENABLE_PYTHON=ON` (in the library section of `CMakeLists.txt`,
also with `ENABLE_PROGRAM=OFF`) into `<build>/python`, not installed (pip/conda packaging: after the validation, the
user's decision). Its test, `ctest -R Python` (`PythonModule`: `library/python/tests/test_module.py`, the toy system of
`UAMGLibraryTest.cpp` against SciPy's direct solver). No NGSolve in the repository (the user's decision, 2026-10-19: no
new dependency): its test is `~/code/fhhos4-ngsolve/test_ngsolve.py`, run by hand.
Environment: conda env `fhhos4-python` (python 3.12, numpy, scipy, nanobind, cmake, ninja, cxx-compiler, eigen 5.0.1;
NGSolve from pip, not on conda-forge), build directory `build-python` (`-DENABLE_PROGRAM=OFF -DENABLE_PYTHON=ON`, ~1 min).

NGSolve (6.2.2608; first target of the module; DOLFINx after work on MPI in fhhos4, the user's decision): its HDG
(`L2 * FacetFESpace`, primal, SPD) fits the uncondensed interface. `test_ngsolve.py` (outside the repository): `a.mat.CSR()`
(views of NGSolve's arrays: uint64 row pointers, int32 columns, read in place), the DoF lists from `GetDofNrs()`, the
Dirichlet faces left out (`fes.FreeDofs()`: all or none of a face's DoFs), `ExtractSymmetricBlocks()`, `Gather()`,
`Scatter()`; `condense=False` (with `condense=True`, `a.mat` keeps only the Schur complement, the cell blocks are not
stored). Solutions bitwise those of the former extraction by SciPy (2026-10-19). NGSolve's sparse Cholesky (the
reference) is not deterministic: its error varies at 1e-14 from run to run.
phi_0 = 1 (interpolation of 1: ones). Bases not ordered by degree (Dubiner's order: 2D cell k=3: degrees 0 1 2 3 1 2 3 2
3 3; 3D cells and faces too, 2D faces ordered; Q^k on quads/hexes, Q-degrees not ordered either): since 2026-10-20 (the
user's decision) the DoF lists sort the DoFs of each cell and face by degree (`ByDegree()`: the degrees measured by
least-squares fits of `GetFE().CalcShape()`, the same order on every element), default p-coarsening: k <= 2 unchanged,
k=3 21 / 69 iterations instead of 44 / 96 with `ToZero` (Square-tri 940 / Cube-tet 1259 cells, tolerance 1e-12), k=4
23 / 59. Solution = NGSolve's sparse Cholesky to 1e-11..1e-14 (2D, 3D, k = l = 1..4; k = 1, l = 2). Iterations with `ToZero`
(tolerance 1e-8, 2026-10-19), U-AMG `bcheb4` / `bgs` / NGSolve's BDDC + CG, penalty alpha/h: Square-tri 3708 cells, alpha = 4k^2
(NGSolve's tutorial): k=1 11/14/12, k=2 14/19/15, k=3 18/22/19; alpha = 10(k+1)^2: 24-29 / 32-36 / 21-30. Cube-tet 3937
cells, alpha = 5(k+1)^2: k=1 24/28/28, k=2 39/44/45, k=3 46/53/55 (alpha = 4 at k=1 is not SPD in 3D: the Chebyshev
setup refuses it). The counts follow the penalty for both solvers.
Times against NGSolve's BDDC + its CG (2026-10-19, `build-python/bench/`, not in the repository; the same configured
tolerance 1e-8, each solver its own norm, the user's decision: BDDC's true residuals 6e-8..1.4e-6, U-AMG's <= 1e-8;
straight to degree 0; k = 3 skipped, the user's decision): BDDC / U-AMG setup + solve on 8 threads, Square-tri
(alpha = 4k^2) 92k-370k cells: k=1 1.27-1.33, k=2 0.80-0.90; Cube-tet (alpha = 5(k+1)^2): k=1 25k cells 1.17, 69k cells
1.96 (BDDC's setup 1.70 s vs 0.43 s), k=2 9k-27k cells 0.67-0.76. 1 thread: the same trends. With the new interface (2026-10-19,
`run_bddc.sh`, `summary_bddc.py`, the same session, best of 2), 8 threads: setup + solve k=1 1.35-1.39 (2D),
1.34-2.22 (3D), k=2 0.96-1.02 (2D), 0.74-0.92 (3D); end to end (U-AMG: uncondensed assembly + extraction + setup + solve;
BDDC: condensed assembly + its setup + CG + cell unknowns) k=1 1.11-1.82, k=2 0.69-0.87 (1 thread: 1.16-1.92,
0.73-0.95). At k=2 the solve is the gap (16 vs 13-14 iterations in 2D, 42-43 vs 49-52 in 3D, but costlier ones). k=3 (sorted by degree, 2026-10-20,
`run_k3.sh`): setup + solve 0.81-0.85 (2D), 0.69-0.70 (3D), end to end 0.64-0.75 (1 thread: 0.61-0.78), with fewer
iterations (10 vs 18 in 2D, 32-35 vs 54-64 in 3D). BDDC with an AMG coarse
solver (`coarsetype="h1amg"`) always slower (51-91 iterations); NGSolve's sparse Cholesky a tie in 2D, 6x slower than
U-AMG in 3D at 25k cells. The extraction of the blocks (2026-10-19, `ExtractSymmetricBlocks()` on NGSolve's arrays
read in place; before: SciPy's slicing and the conversion to Eigen): U-AMG's total -5 to -26%, `PERFORMANCE.md`.

### U-AMG: details and measurements (kept out of `README_library.md`)

- Condensation by cells (setup from the blocks): computes the rows of `A` in the pattern of `A_FF`, completed (and
  sorted) only in the rows where `A_FF` lacks couplings of the faces of a cell (zeros dropped by an assembly): 1.3x (1
  thread) to 1.8x (8 threads) slower when all the rows are completed. HArDCore3D 32^3, k = 1 (285,696 face DoFs, 9.2 M
  non-zeros in `A`): setup from the blocks 0.71-0.80 s vs 0.38-0.50 s from `A` on 1 thread, 0.25-0.35 s vs 0.17-0.29 s
  on 8; same solve. The setups from the blocks keep the inverses of the blocks of `A_TT` and `A_TT^{-1} A_TF`, and
  read `A_FT` at each solve: the transpose of `A_TF` (symmetric), or the caller's in place (non-symmetric).
- p-levels: `P_CS` Minus2 (k, k-2, ..., 0), blocks extracted (first coefficients of each block,
  `SetupOperatorByBlockExtraction()`), prolongation = injection of the low modes. `PCoarsening::ToZero` (`-p-cs =0`, k then 0,
  for bases not ordered by degree, e.g. NGSolve's): more iterations from k = 3 (FCG, 2026-10-19; `bgs,rbgs` / `bcheb4`):
  Square-tri n=64 k=3 10 -> 19 / 14 -> 27, k=4 9 -> 21 / 14 -> 30, k=5 7 -> 23 / 11 -> 32; Cube-tet n=8 k=3 29 -> 37 /
  50 -> 60. h-levels: multiple pairwise
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
- Eigen and the interface (2026-10-19, replaces the layout check of 2026-10-17): no Eigen object crosses the
  boundary of the library, any version of Eigen in the client. `BorrowedSparseMatrix` relies on the protected members of
  Eigen's `SparseMatrix` and `CompressedStorage`, the same in Eigen 3.4.0 and 5.0.1 (`tests/BorrowedSparseMatrixTest.cpp`;
  checked under valgrind with both). Eigen 3.4 only lacked `saveMarketDense` (`Utils/ExportModule.h` writes the format
  itself). Eigen 5's legacy macros say 3.5.0 (`EIGEN_WORLD_VERSION` stays 3; `EIGEN_PATCH_VERSION` exists from 5.0).
- Client codes with column-major matrices (2026-10-19, the user: minimal changes in the clients): HArDCore3D's bindings
  build their blocks row-major (the same parallel builders, by rows), its symmetric `A` through `TransposedView()`, its
  non-symmetric `SysMat` through `ToCSR()`; Disk++: `RowMajor` on its `ATT`/`ATF`, `TransposedView(assm.LHS)`. Same
  iterations; bitwise results with `ToCSR()`, the errors at 1e-10 relative with `TransposedView()` (its `A^T`).

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
library alone with 64-bit indices (`SMALL_INDEX=OFF`) and the conda test client, and in a parallel job the conda
package with its tests, then creates the GitHub release "Release X.Y.Z" (notes: the message of the tag). It pins GMSH,
Eigen, CGAL, the compiler and rattler-build to the versions the tests were validated with: update them with the local
environment (envs `fhhos4` and `rattler`). After a change of `release.yml` or of its pins, run it by hand on master
(Actions > Release > Run workflow: no release) before pushing the tag; otherwise push the tag directly. After the release, conda-forge's bot updates the feedstock; the Spack package
needs a PR to spack-packages (`spack checksum --add-to-package fhhos4 X.Y.Z`, then `spack install --test=root` and
`spack test run`).

## Documentation (README files)

Three READMEs, each for its readers; keep them up to date in the same commit as the change they describe:
- `README.md`: minimal. What fhhos4 is, its two parts (one paragraph each, linking to the two files below), and the
  license. Nothing else.
- `README_program.md`: the program, for its users (researchers running experiments): dependencies, conda environment,
  build (options table: every `option()` of `CMakeLists.txt`), installation, where the program reads and writes,
  troubleshooting.
- `README_library.md`: the library, for the developers of client codes, who don't need the program or the sources. Its
  numbered sections (the user's structure, 2026-10-18): 1. Installation (install from source, the library alone; the
  conda package; link, CMake package), 2. The algorithm in a nutshell (in a few lines: p-levels, then h-levels and how
  the blocks give the aggregation; applicability), 3. How to use (the two interfaces, "Interface 1" the uncondensed
  system, recommended, and "Interface 2" the condensed one, with examples; non-symmetric; Python; inputs, parameters,
  constraints, license, version last). Rules:
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
