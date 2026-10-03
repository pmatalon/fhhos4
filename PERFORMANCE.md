# Performance

Performance is key for this project. This file records what has been measured, how to measure and validate a change,
and the next optimization ideas with their expected gains, so that each optimization starts from data.
Update it after measuring or optimizing something.

## Rules

- **The papers stay reproducible.** With the default options, the iteration counts of `reproducibility/*.md` must not
  change (`ctest` checks them). Aim for changes that are bit-identical: same matrices, same residuals. A change of
  numerics (e.g. another smoother) must be an opt-in option, never the default.
- **Measure before and after**, sequential (`-threads 1`) and parallel (default: all logical cores), elapsed time,
  best of 2 runs. The noise is about ±10% (WSL2). Run nothing else meanwhile: no build, one heavy job at a time
  (13 GB of RAM).
- The 2023 AMG paper reports **sequential** CPU times (§4.2, Fig. 4.5): parallel gains don't show there,
  sequential ones do.

## Machine

AMD Ryzen 7 7735HS (8 cores, 16 threads), 13 GB of RAM, WSL2. No `perf` in WSL2; `gdb` and `valgrind` are not
installed but may be.

## How to measure

The scripts and headers are in `scripts/perf/` (run the scripts from `build/`, with the conda env activated; keep a
copy of the old binary, e.g. `cp bin/fhhos4 /tmp/fhhos4_before`, to compare with it).

- **Times**: `scripts/perf/timing.sh <output dir> [binary]` runs the reference cases below, sequential (`-threads 1`)
  and parallel (`-threads 0`), twice each, and writes the best times in `best.txt`. They come from the `Setup` and
  `Solving` rows (`Elapsed time` column) of the table printed after the solve.
- **Reference cases** (U-AMG):
  - Cube-tet, k=0: `./bin/fhhos4 -geo cube -mesh tetra -k 0 -n 32 -s fcguamg -no-cache` (318k face unknowns)
  - Cube-cart-aniso100: `./bin/fhhos4 -geo cube -mesh cart -mesher inhouse -k 0 -n 64 -aniso 100 -s fcguamg -no-cache`
  - High order: `./bin/fhhos4 -geo cube -mesh tetra -k 2 -n 16 -s fcguamg -no-cache` (262k face unknowns)
  - The variants of Fig. 4.5 / Table 4.4: `-coarsening-prolong 6` (P_F, default), `4` (Q_F^smooth), `5` (P_F^(0)), `3` (Q_F).
- **Profiling**: no sampling profiler, so temporary wall-clock timers around the steps, per multigrid level:
  `scripts/perf/Prof.h` (`PROF_START`/`PROF_STOP`, printed at exit). Include it in the files to instrument, never
  commit the instrumentation.
- **Validation of a change**:
  1. `ctest` (85 tests).
  2. `scripts/perf/amg_papers.sh <output dir> [binary]`: the AMG configurations of the papers, at reduced sizes (the
     paper sizes don't fit in 13 GB), 30 runs in ~6 minutes. Then `scripts/perf/fingerprint.sh <dir before> <dir after>`
     compares their iteration tables (iteration, inner iterations, residual): they must be identical.
     The runs: Cube-cart n=32, Cube-tet n=16, Complex-tet (`-geo platewith4holes`) n=4, Heterog1e8 n=256,
     Cube-cart-aniso100 n=32, Cube-tet-aniso20 n=16 and Cube-tet n=32 (Fig. 4.4), each with `fcguamg` and
     `fcgaggregamg`; Cube-tet n=16 with `-cs mpa|dpa`; Cube-tet n=16 and Cube-cart-aniso100 n=32 with
     `-coarsening-prolong 3|4|5|6`; the biharmonic Table 3 runs (n=8,16, `-bihar-prec s|no`) and the commented-out
     k=3 runs (`-opt2 0|2`) of `2023_biharmonic_problem.md`.
  3. When replacing a computation: temporary bit-for-bit checks of each new matrix against the old computation,
     with `scripts/perf/BitCheck.h` (same structure, `memcmp` of the values), on the runs of `amg_papers.sh`
     (`grep BITCHECK` in the logs).

## Exactness facts (to keep the results bit-identical)

- **Eigen's sparse product** `A * B` (row-major): each coefficient is summed in the storage order of the row of `A`,
  the first term is assigned (not added to 0), the explicit zeros are kept, and the rows are sorted.
  `SparseMatrixOps::Multiply` does the same, whatever the number of threads. Eigen evaluates the products from left
  to right: `P.transpose() * S * P` is `Multiply(Multiply(Pt, S), P)`.
- **`S.selfadjointView<Eigen::Lower>()`** only reads the lower triangle: `SparseMatrixOps::FullFromLower(S)`. It
  matters: the Galerkin operators are not bitwise symmetric.
- **`NonZeroCoefficients::Add`** drops the coefficients of absolute value <= `NonZeroCoefficients::ZeroThreshold`
  (1e-15). `setFromTriplets` sorts the rows and sums the duplicates.
- **The hybrid smoothers** (`hbgs`, `hrbgs`: Gauss-Seidel inside the chunk of rows of each thread, Jacobi between
  the chunks) give results that depend on the number of threads.

## Done

- **U-AMG setup** (commit `fe213081`, 2026-10-02): coarse `A_F_F` and chained `Q_F` only computed when the options
  use them; `BlockJacobi::IterationMatrix()` built in parallel, and only for the rows of the removed faces in the
  prolongations 4 and 6; all the sparse products of the setup computed in parallel (`src/Utils/SparseMatrixOps.h`).
  Bit-identical. Setup time, before -> after:

  | Case | Sequential | Parallel (16 threads) |
  |---|---|---|
  | Cube-tet k=0 n=32, P_F | 6.41 -> 3.82 s | 5.66 -> 2.06 s |
  | Cube-tet k=0 n=32, Q_F^smooth / P_F^(0) / Q_F | 6.12 / 4.67 / 4.35 -> 3.53 / 3.47 / 3.30 s | 5.74 / 3.97 / 3.78 -> 1.97 / 2.04 / 1.89 s |
  | Cube-cart-aniso100 n=64 | 8.13 -> 4.78 s | 7.25 -> 2.92 s |
  | Cube-tet k=2 n=16 | 17.79 -> 9.07 s | 17.03 -> 3.91 s |

## Current profile (commit `fe213081`, U-AMG, default OpenMP settings)

| Case | Setup, sequential / parallel | Solve, sequential / parallel |
|---|---|---|
| Cube-tet k=0 n=32 | 3.8 / 2.1 s | 1.0 / 1.0 s |
| Cube-cart-aniso100 n=64 | 4.8 / 2.9 s | 0.6 / 0.7 s |
| Cube-tet k=2 n=16 | 9.1 / 3.9 s | 6.3 / 7.0-7.3 s |

- **The solve doesn't benefit from the threads.** The smoothing is 75-80% of it (0.75 s at k=0, 5.6 s at k=2), and
  the default smoother (lexicographic Gauss-Seidel, as in the paper) is sequential. At k=2, 16 threads are even
  slower than 1, apparently because the idle threads spin on the same cores while the smoother runs (see idea 1).
- **Setup, k=0, 16 threads (2.1 s)**: hybrid meshes 0.9 s (build 0.16 s + auxiliary coarse mesh 0.12 s + coarsening
  0.45 s, including the sequential pairwise aggregation 0.18 s and interface collapsing 0.19 s + freeing 0.17 s);
  sparse products 0.4 s; `FullFromLower(S)` 0.16 s; `BuildQ_F` 0.10 s; `Theta()` 0.08 s; other transposes 0.09 s.
- **Setup, k=2, 16 threads (3.7 s)**: sequential transposes 1.2 s (`FullFromLower(S)` 0.81 s, `R = P^T` 0.24 s,
  column-major copy of `A_T_F` in `HybridAlgebraicMesh::Build()` 0.20 s); products 1.2 s (including `P^T*S` 0.52 s
  and `(P^T*S)*P` 0.43 s); `Theta()` 0.28 s; P chaining 0.20 s; hybrid meshes 0.4 s.

## Next ideas, by expected gain

1. **`OMP_WAIT_POLICY=passive`** — measured, no code, results unchanged. k=2 (16 threads): solve 7.0-7.3 -> 6.2-6.3 s
   (-12%), setup unchanged or slightly better; k=0: no notable change. (`OMP_WAIT_POLICY=active` is much worse:
   solve 8.9 s.) libgomp reads it at startup, so it cannot be set from `main()`: set it in the environment, e.g.
   `conda env config vars set OMP_WAIT_POLICY=passive -n fhhos4`, and document it.
2. **Parallel transpose** — estimated, small effort, bit-identical (a transpose is exact). Replaces Eigen's sequential
   transposes in `FullFromLower`, `UncondensedLevel::SetupRestriction` (`R = P^T`), `HybridAlgebraicMesh::Build()`
   (column-major copy of `A_T_F`) and `SparseMatrixOps::Transpose`. They take 1.2 s of the 3.7 s of setup at k=2 and
   0.25 s of the 2.1 s at k=0: if they scale like the product (5-7x on 8 cores), setup -30% at k=2, -10% at k=0.
   `FullFromLower` only needs the transpose of the strictly lower part.
3. **Flat data structures for the hybrid algebraic meshes** (`HybridAlgebraicMesh`) — estimated, medium-large
   effort, main lever at k=0 (the paper's case). Vectors of vectors, a mutex per element, a `std::map` of neighbours
   per aggregate, copied vectors: building, coarsening and freeing the meshes take 0.9 s of the 2.1 s of setup at
   k=0. CSR-like arrays and prefix sums (same numbering of the aggregates and coarse faces, so bit-identical) might
   halve it: setup ~2.1 -> ~1.6 s. The pairwise aggregation itself (greedy, priority order) must stay sequential:
   parallelizing it changes the aggregates.
4. **Parallel smoother, as an option** — measured with the existing hybrid block Gauss-Seidel
   (`-smoothers hbgs,hrbgs`, 16 threads, passive wait policy): k=2: 40 iterations instead of 35, solve 6.2 -> 3.3 s
   (smoothing 4.8 -> 1.5 s); k=0: 31 iterations instead of 28, solve 0.94 -> 1.40 s (slower: no gain on scalar
   rows). Changes the iteration counts and depends on the number of threads: opt-in only. At k=0, a scalar parallel
   Gauss-Seidel (e.g. multicolour) would be needed; its gain is unknown.

Other candidates, not measured yet:

- `Level::ComputeGalerkinOperator()` computes `R * A * P` with Eigen's sequential products: used by the multigrids
  with the Galerkin operator that don't set their coarse operators themselves (e.g. the geometric multigrid with
  `-g 1`). `SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(R, A), P)` gives the same matrix, bit for bit.
  Profile those setups first.
- C-AMG (`AggregAMG::GalerkinOperator`): sequential loop and triplets. Making it faster changes the U-AMG/C-AMG
  comparison of the paper (whose timings are tied to release 1.0 anyway).
- `HybridAlgebraicMesh::Theta()` and `CouplingValue()` extract dense blocks with `SparseMatrix::block()`, which scans
  the row from its start: the cost grows with the block size, i.e. with k.
- `BlockJacobi::Setup()` factorizes all the diagonal blocks, while the prolongations 4 and 6 only need those of the
  removed faces (0.02 s at k=0, 0.09 s sequential at k=2).
- Memory: peak 1.33 GB for Cube-tet k=0 n=32 (whole run). The paper's n=64 extrapolates to 10-11 GB, too close to
  13 GB. Running the paper sizes here would first need a memory profile (probably the mesh and the assembly rather
  than the AMG: unverified).
