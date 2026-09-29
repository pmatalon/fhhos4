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

Tests don't go through the CLI as a subprocess. Instead, `tests/support/RunHelper.h` builds a `ProgramArguments` struct directly, applies the same argument-defaulting cascade as `main.cpp` (`ApplyProgramArgumentDefaults`), and calls `Program_Diffusion_HHO<Dim>::Execute()` in-process. So every test case below is annotated with its CLI-equivalent `fhhos4 ...` command for readability, even though the test itself never spawns `fhhos4`.

There is no golden/reference-file comparison anywhere in the suite — every assertion is computed in-code: a least-squares convergence-order slope, a generous iteration-count bound, or (for two intentionally-broken configurations) a death-test exit code.

## How to run

```bash
cd build
ctest                                      # run the whole suite
ctest -R Kellogg                           # run a subset by name (regex)
./bin/fhhos4_tests --gtest_filter=*Kellogg*  # or filter the binary directly
```

The full suite is 39 CTest cases and runs in about 2 minutes serially (the slowest single case, `MeshIndependenceTest/("stri", 3)`, takes ~19s).

## `HMultigrid2020Test.cpp`

Validates *"An h-multigrid method for Hybrid High-Order discretizations"* (Di Pietro, Hulsemann, Matalon, Mycek, Rude, Ruiz, SIAM J. Sci. Comput. 2021) — see `reproducibility/2020_MG_for_HHO.md`.

- **`Square/ConvergenceOrderTest.MatchesTheoreticalOrder`** (8 cases: `{cart, stri} × k=0..3`) — the L2 error should converge at the theoretical rate (h² for k=0, h^(k+2) for k≥1) as the mesh is refined n=8,16,32,64.
  CLI: `fhhos4 -geo square -mesh {cart|stri} -mesher inhouse -s mg -cycle V,1,1 -k {0..3} -n {8,16,32,64}`
  MG iteration counts observed per n=8/16/32/64 (n=8/16 sit in a degenerate pre-asymptotic regime where the multigrid hierarchy barely has any levels, hence the 1's):
  - cart: k=0 → 1,1,6,7 · k=1 → 1,1,8,8 · k=2 → 1,9,10,10 · k=3 → 1,10,11,11
  - stri: k=0 → 1,1,10,10 · k=1 → 1,11,12,12 · k=2 → 1,10,10,10 · k=3 → 1,10,10,10

- **`Square/MeshIndependenceTest.IterationCountStaysBounded`** (8 cases: `{cart, stri} × k=0..3`) — the V(1,1)-cycle iteration count should stay essentially mesh-independent as n=32→64→128. Assertion: largest ≤ 2×smallest and largest ≤ 40 (a generous bound rather than exact constancy).
  Observed: cart: k=0 → 6,7,7 · k=1 → 8,8,8 · k=2 → 10,10,10 · k=3 → 11,11,11; stri: k=0 → 10,10,10 · k=1 → 12,12,12 · k=2 → 10,10,9 · k=3 → 10,10,10.
  Note: the in-code comment justifying the generous bound cites an older observation of "19 → 24 → 27 for the unstructured mesh at k=0", but the current numbers above (stri k=0: 10,10,10) are flat — the bound still holds comfortably either way.

- **`HMultigrid2020.DegradesAsDocumented_StructuredTetra_K0`** (1 case) — the paper's Figure 4.1 documents this configuration (`cube`/`stetra`, k=0, V(2,2)) as diverging: with the standard coarsening, the convergence rate degrades toward 1. Assertion: more than 60 iterations (observed: 122 at n=16, vs. 24 for k=1).
  CLI: `fhhos4 -geo cube -mesh stetra -mesher inhouse -s mg -cycle V,2,2 -k 0 -n 16`

- **`SquareFourQuadrants/KelloggTest.ConvergesWithinBound`** (4 cases: k=0..3) — the heterogeneous Kellogg benchmark (Figure 4.7) should converge within 30 iterations at n=256.
  CLI: `fhhos4 -geo square4quadrants -tc kellogg -mesh cart -mesher inhouse -s mg -cycle V,1,1 -k {0..3} -n 256`
  Observed iteration counts: k=0 → 10 · k=1 → 9 · k=2 → 10 · k=3 → 10.

- **`SquareFourQuadrants/HeterogeneityRatioSweepTest.IterationCountRatioBounded`** (4 cases: k=0..3) — iteration count should stay bounded as the `-heterog` diffusion-coefficient ratio sweeps 1e0/1e2/1e4/1e6/1e8, with the Galerkin operator enabled (Figure 4.9(a)). n=64, V(0,3). Assertion: max ≤ 3×min.
  Observed: k=0 → 7,7,7,7,7 · k=1 → 9,9,9,9,9 · k=2 → 10,10,10,10,10 · k=3 → 11,11,11,11,11 — completely flat at these settings; the ratio has essentially no effect.

## `HpStrategies2022Test.cpp`

Validates *"High-order multigrid strategies for HHO discretizations of elliptic equations"* (Di Pietro, Matalon, Mycek, Rude, Numer. Linear Algebra Appl. 2022) — see `reproducibility/2022_high_order_strategies.md`.

- **`HpStrategies2022.BasisNormalization_OrthonormalDiverges`** / **`BasisNormalization_OrthogonalConverges`** (2 cases) — §3.4.1: with local refinement, orthonormalized element bases (`-e-ogb 3`) make the multigrid diverge (death test, `EXIT_FAILURE`), while orthogonal bases (`-e-ogb 1`) converge (observed: 13 iterations; assertion ≤ 30).
  CLI: `fhhos4 -geo square4quadrants_tri_localref -no-cache -tc square -cs r -k 1 -n 32 -e-ogb {3|1}`

- **`SquareCart/HPConfigTest.AllConverge`** (8 cases: `{mg, fcgmg} × hp-config{1,2,3,4}`) — every hp-multigrid coarsening strategy (h-only, p→h, p→h with h-prolongation, hp→h) should converge at high order (k=5, n=32), both as a stand-alone solver and as an FCG preconditioner. Assertion: 0 < iterations ≤ 100.
  Observed iteration counts: `mg` → hp1=15, hp2=7, hp3=10, hp4=15; `fcgmg` → hp1=11, hp2=6, hp3=9, hp4=11. hp-config 2 (p→h, injection/remove-higher-orders) is consistently the fastest of the four.

- **`SquareCart/ConvergenceOrderHighOrderTest.MatchesTheoreticalOrder`** (4 cases: k=2..5) — at high order, with hp-config 2 and a tight tolerance (1e-12), the L2 error should still follow the theoretical h^(k+2) order as n=16→32.
  Observed MG iteration counts: k=2 → 5,7 · k=3 → 9,9 · k=4 → 7,7 · k=5 → 7,7.

## Coverage summary

- 39 CTest cases total across the 2 files above, all exercising `Program_Diffusion_HHO` — i.e. **diffusion (Poisson-type) problems, HHO discretization, static condensation**.
- Dimensions: 2D (37 cases) plus one 3D death test that doesn't even reach the solver. No 1D coverage (and `RunDiffusionHHO` throws for it, consistent with `ENABLE_1D=OFF` by default).
- Meshes: in-house `cart`, `stri`, `stetra` and GMSH `cart` (the 2022 tests leave the mesher at its default, which resolves to GMSH). No `poly`/agglomerated (CGAL) coverage.
- Solvers: `mg` and `fcgmg` only. No coverage of `lu`, `ch`, `cg`, `eigencg`, `uamg`, `aggregamg`, `agmg`, or `p_mg`.
- Test cases (`-tc`): the default `sine` test case on `square`, `kellogg` on `square4quadrants`, a homogeneous `-heterog` ratio sweep on `square4quadrants`, and one locally-refined-mesh case.
- **Not covered**: biharmonic (`-pb bihar` / `bihardd`), DG/FEM discretizations, full-Neumann BCs, anisotropy, non-condensed systems, and most solver/mesh codes listed above.
- Two tests are *intentionally* death tests that regression-document known, currently-failing configurations drawn straight from the papers' own reproduction recipes — they are not bugs to silently "fix"; changing their behavior should be a deliberate, separate change.
- Papers with reproducibility docs but no corresponding tests yet: `2021_non_nested_MG_for_HHO.md`, `2023_AMG_for_hybrid_methods.md`, `2023_biharmonic_problem.md`, `2024_HHO_HDG_demo_framework.md`.

## Adding a test

Drop a new `*.cpp` file under `tests/` — it's picked up automatically by the CMake glob, no build-file changes needed. Reuse the helpers in `tests/support/RunHelper.h`:
- `fhhos4_tests::RunDiffusionHHO(ProgramArguments)` — applies CLI-equivalent defaulting and runs a diffusion HHO problem in-process, returning iteration count and L2 error.
- `fhhos4_tests::SyncGlobalProgramState<Dim>(args)` — synchronizes global state (`Utils::ProgramArgs`, mesh directories, GMSH cache flags) that mesh-construction code reads directly; called automatically by `RunDiffusionHHO`.
- `fhhos4_tests::EstimateConvergenceOrder(h, errors)` — least-squares log-log slope, for convergence-order assertions.

A run that's expected to diverge (`Utils::FatalError`) terminates the process, so wrap the call in GTest's `EXPECT_EXIT`/`ASSERT_EXIT` rather than a plain assertion.
