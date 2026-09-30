# tests

This directory contains the validation test suite for `fhhos4`, built as a single GoogleTest binary (`fhhos4_tests`). It is wired up at the end of the top-level `CMakeLists.txt`:

```cmake
enable_testing()
find_package(GTest REQUIRED)
include(GoogleTest)

file(GLOB_RECURSE test_source_files CONFIGURE_DEPENDS tests/*.cpp)
add_executable(fhhos4_tests ${test_source_files})
target_include_directories(fhhos4_tests PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(fhhos4_tests PRIVATE fhhos4_core GTest::gtest_main)
gtest_discover_tests(fhhos4_tests)
```

Any new `tests/*.cpp` file is picked up automatically (`CONFIGURE_DEPENDS`) — no CMake edit needed to add a test.

Tests don't go through the CLI as a subprocess. Instead, `tests/support/RunHelper.h` builds a `ProgramArguments` struct directly, applies the same argument-defaulting cascade as `main.cpp` (`ApplyProgramArgumentDefaults`), and calls `Program_Diffusion_HHO<Dim>::Execute()` in-process. So every test case (here, and in a comment above each test in the source) is annotated with its CLI-equivalent `./bin/fhhos4 ...` command, to run it manually from `build/`, even though the test itself never spawns `fhhos4`.

Assertions are of three kinds: a least-squares convergence-order slope; an iteration count, either exact (the 2020 paper's published values) or bounded (the 2022 paper); or (for two intentionally-broken configurations) a degradation or a death-test exit code.

## How to run

```bash
cd build
ctest                                      # run the whole suite
ctest -R Kellogg                           # run a subset by name (regex)
./bin/fhhos4_tests --gtest_filter=*Kellogg*  # or filter the binary directly
```

The full suite is 31 CTest cases and runs in under 30 seconds serially (no single case takes more than ~2s). Tests are deliberately kept small (N ≤ 64, k ≤ 1 for the multigrid tests).

## `HHOConvergenceTest.cpp`

A-priori convergence of the HHO discretization, independent of any solver study: the linear system is solved directly (sparse Cholesky, `-s ch`), so the L2 error is the discretization error.

- **`Square/ConvergenceOrderTest.MatchesTheoreticalOrder`** (8 cases: `{cart, stri} × k=0..3`) — the L2 error should converge at the theoretical rate (h² for k=0, h^(k+2) for k≥1) as the mesh is refined n=8,16,32. Assertion: least-squares slope within ±0.3.
  CLI: `./bin/fhhos4 -geo square -mesh {cart|stri} -mesher inhouse -s ch -k {0..3} -n {8|16|32}`

## `HMultigrid2020Test.cpp`

Validates *"An h-multigrid method for Hybrid High-Order discretizations"* (Di Pietro, Hulsemann, Matalon, Mycek, Rude, Ruiz, SIAM J. Sci. Comput. 2021) — see `reproducibility/2020_MG_for_HHO.md`.

The expected iteration counts are **exactly** the paper's, taken from the CSV files of the accepted version (`Multigrid_for_HHO_SISC.zip`, `Results/`; in the file names, `p` = k+1). Only the smallest mesh sizes and k=0,1 are tested, to keep the suite fast.

- **`Square/IterationCountTest.MatchesPaper`** (4 cases: `{cart, stri} × k=0,1`) — Figure 4.1, V(1,1) cycle, n=32,64. Reference: `2D_scalability_homogeneous_V11_g0_p{1|2}_{cart|tri}.csv`.
  CLI: `./bin/fhhos4 -geo square -mesh {cart|stri} -mesher inhouse -s mg -cycle V,1,1 -k {0|1} -n {32|64}`
  Expected: cart: k=0 → 13,14 · k=1 → 16,18; stri: k=0 → 19,24 · k=1 → 24,24.
  All other points of the CSVs checked up to n=128 (cart/stri, k=0..3) match, except stri, k=1, n=128 (not tested): the code gives 25 iterations, the CSV reports 24. The CSV value is most likely wrong: every code version from the introduction of standard coarsening for triangles (Oct 2019) through the paper's submission (the CSVs are unchanged since the June 2020 submission) up to today gives 25, with the exact same residual history (1.44e-8 at iteration 24, far from the 1e-8 tolerance). The rest of the curve matches (n=256, 512 → 26, 26).

- **`HMultigrid2020.DegradesAsDocumented_StructuredTetra_K0`** (1 case) — Figure 4.1 documents this configuration (`cube`/`stetra`, k=0, V(2,2)) as diverging (its CSV has no iteration count): with the standard coarsening, the convergence rate degrades toward 1. Assertion: more than 60 iterations (observed: 122 at n=16, vs. 24 for k=1).
  CLI: `./bin/fhhos4 -geo cube -mesh stetra -mesher inhouse -s mg -cycle V,2,2 -k 0 -n 16`

- **`SquareFourQuadrants/KelloggTest.MatchesPaper`** (2 cases: k=0,1) — Figure 4.7, heterogeneous Kellogg benchmark, V(1,1) cycle, n=32,64. Reference: `Kellogg_scalability_V11_g0_p{1|2}_cart.csv`.
  CLI: `./bin/fhhos4 -geo square4quadrants -tc kellogg -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k {0|1} -n {32|64}`
  Expected: k=0 → 13,13 · k=1 → 15,17.

- **`SquareFourQuadrants/HeterogeneityRatioSweepTest.MatchesPaper`** (2 cases: k=0,1) — Figure 4.9(a): with the Galerkin operator and the heterogeneous weighting, the V(0,3) iteration count doesn't depend on the `-heterog` ratio (1e0/1e2/1e4/1e6/1e8), n=64. Reference: `2D_heterogeneity_chiasmus_n64_V03_g1_cart.csv`.
  CLI: `./bin/fhhos4 -geo square4quadrants -mesh cart -mesher inhouse -n 64 -s mg -g 1 -cycle V,0,3 -k {0|1} -heterog {1e0|1e2|1e4|1e6|1e8}`
  Expected: k=0 → 7 · k=1 → 9, for every ratio.

## `HpStrategies2022Test.cpp`

Validates *"High-order multigrid strategies for HHO discretizations of elliptic equations"* (Di Pietro, Matalon, Mycek, Rude, Numer. Linear Algebra Appl. 2022) — see `reproducibility/2022_high_order_strategies.md`.

- **`HpStrategies2022.BasisNormalization_OrthonormalDiverges`** / **`BasisNormalization_OrthogonalConverges`** (2 cases) — §3.4.1: with local refinement, orthonormalized element bases (`-e-ogb 3`) make the multigrid diverge (death test, `EXIT_FAILURE`), while orthogonal bases (`-e-ogb 1`) converge (observed: 13 iterations; assertion ≤ 30).
  CLI: `./bin/fhhos4 -geo square4quadrants_tri_localref -no-cache -tc square -cs r -k 1 -n 32 -e-ogb {3|1}`

- **`SquareCart/HPConfigTest.AllConverge`** (8 cases: `{mg, fcgmg} × hp-config{1,2,3,4}`) — every hp-multigrid coarsening strategy (h-only, p→h, p→h with h-prolongation, hp→h) should converge at high order (k=5, n=32), both as a stand-alone solver and as an FCG preconditioner. Assertion: 0 < iterations ≤ 100.
  CLI: `./bin/fhhos4 -geo square -mesh cart -cs r -k 5 -n 32 -s {mg|fcgmg} -tol 1e-10 -hp-config {1|2|3|4}`
  Observed iteration counts: `mg` → hp1=15, hp2=7, hp3=10, hp4=15; `fcgmg` → hp1=11, hp2=6, hp3=9, hp4=11. hp-config 2 (p→h, injection/remove-higher-orders) is consistently the fastest of the four.

- **`SquareCart/ConvergenceOrderHighOrderTest.MatchesTheoreticalOrder`** (4 cases: k=2..5) — at high order, with hp-config 2 and a tight tolerance (1e-12), the L2 error should still follow the theoretical h^(k+2) order as n=16→32.
  CLI: `./bin/fhhos4 -geo square -mesh cart -cs r -k {2..5} -n {16|32} -s fcgmg -tol 1e-12 -hp-config 2`
  Observed MG iteration counts: k=2 → 5,7 · k=3 → 9,9 · k=4 → 7,7 · k=5 → 7,7.

## Coverage summary

- 31 CTest cases total across the 3 files above, all exercising `Program_Diffusion_HHO` — i.e. **diffusion (Poisson-type) problems, HHO discretization, static condensation**.
- Dimensions: 2D (30 cases) plus one 3D case (`stetra`, k=0). No 1D coverage (and `RunDiffusionHHO` throws for it, consistent with `ENABLE_1D=OFF` by default).
- Meshes: in-house `cart`, `stri`, `stetra` and GMSH `cart` (the 2022 tests leave the mesher at its default, which resolves to GMSH). No `poly`/agglomerated (CGAL) coverage.
- Solvers: `mg`, `fcgmg` and `ch` only. No coverage of `lu`, `cg`, `eigencg`, `uamg`, `aggregamg`, `agmg`, or `p_mg`.
- Test cases (`-tc`): the default `sine` test case on `square`, `kellogg` on `square4quadrants`, a homogeneous `-heterog` ratio sweep on `square4quadrants`, and one locally-refined-mesh case.
- **Not covered**: biharmonic (`-pb bihar` / `bihardd`), DG/FEM discretizations, full-Neumann BCs, anisotropy, non-condensed systems, and most solver/mesh codes listed above.
- Two tests *intentionally* document known, failing configurations (a death test in the 2022 file, a near-divergence in the 2020 file) drawn straight from the papers' own reproduction recipes — they are not bugs to silently "fix"; changing their behavior should be a deliberate, separate change.
- Papers with reproducibility docs but no corresponding tests yet: `2021_non_nested_MG_for_HHO.md`, `2023_AMG_for_hybrid_methods.md`, `2023_biharmonic_problem.md`, `2024_HHO_HDG_demo_framework.md`.

## Adding a test

Drop a new `*.cpp` file under `tests/` — it's picked up automatically by the CMake glob, no build-file changes needed. Reuse the helpers in `tests/support/RunHelper.h`:
- `fhhos4_tests::RunDiffusionHHO(ProgramArguments)` — applies CLI-equivalent defaulting and runs a diffusion HHO problem in-process, returning iteration count and L2 error.
- `fhhos4_tests::SyncGlobalProgramState<Dim>(args)` — synchronizes global state (`Utils::ProgramArgs`, mesh directories, GMSH cache flags) that mesh-construction code reads directly; called automatically by `RunDiffusionHHO`.
- `fhhos4_tests::EstimateConvergenceOrder(h, errors)` — least-squares log-log slope, for convergence-order assertions.

A run that's expected to diverge (`Utils::FatalError`) terminates the process, so wrap the call in GTest's `EXPECT_EXIT`/`ASSERT_EXIT` rather than a plain assertion.
