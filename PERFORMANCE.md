# Performance

Performance is key for this project. This file records what has been measured, how to measure and validate a change,
and the next optimization ideas with their expected gains, so that each optimization starts from data.
Update it after measuring or optimizing something.

## Rules

- **The papers stay reproducible.** With the default options, the iteration counts of `reproducibility/*.md` must not
  change (`ctest` checks them).
- **Rounding-level differences are allowed.** A change may reorder the arithmetic operations (the order of a sum,
  several partial sums, vectorized reductions...) if it helps performance, even if the results are then no longer
  bit-identical. A change of the numerical method (another smoother, prolongation, precision...) must be an opt-in
  option, never the default. A bit-identical change is the easiest to validate (`solutions.sh`): prefer it when it
  costs nothing. The sign of a zero never counts (+0 and -0 are equal).
- **Measure before and after**, sequential (`-threads 1`) and parallel (default: all logical cores, 16 here, and since
  2026-10-08 one thread per physical core, 8 here, in the solver; the solver measurements before that date used 16),
  elapsed time,
  best of 2 runs. The noise is about ±10% (WSL2). Run nothing else meanwhile: no build, one heavy job at a time
  (13 GB of RAM).
- The 2023 AMG paper reports **sequential** CPU times (§4.2, Fig. 4.5): parallel gains don't show there,
  sequential ones do.

## Machine

AMD Ryzen 7 7735HS (8 cores, 16 threads), 13 GB of RAM, WSL2. No `perf` in WSL2; `valgrind` is installed (system,
`/usr/bin/valgrind`, not in the conda env), `gdb` is not but may be.

Memory: WSL can crash (the whole VM restarts). It crashed twice on 2026-10-04 during a measurement, cause not
identified: once during runs of `timing.sh`, once during a build with two compilations in parallel (`Program.cpp`
alone peaks at ~3.5-4 GB). The measurement was then redone with `ninja -j1` and a watchdog that kills the compiler
or the run if `MemAvailable` drops below 2.5 GB: the minimum was 7.4 GB (runs of `timing.sh`: at most 2.5 GB each),
no crash. `/tmp` is a tmpfs (in RAM): put no build tree there.

Fresh memory is slow (page faults): the first write to newly allocated memory runs at ~1.3 GB/s, against ~12 GB/s
once mapped, and ~2.5 GB/s with `madvise(MADV_HUGEPAGE)` (transparent huge pages are in `madvise` mode). glibc maps
the large allocations afresh each time (and unmaps them when freed): a large array allocated at each step (e.g. each
coarsening pass) pays it each time. Reuse such buffers (`LocalMatrices` does).

## How to measure

The scripts and headers are in `scripts/perf/` (run the scripts from `build/`, with the conda env activated; keep a
copy of the old binary, e.g. `cp bin/fhhos4 /tmp/fhhos4_before`, to compare with it).

- **Times**: `scripts/perf/timing.sh <output dir> [binary]` runs the reference cases below, sequential (`-threads 1`)
  and parallel (`-threads 0`), twice each, and writes the best times in `best.txt`. They come from the `Setup` and
  `Solving` rows (`Elapsed time` column) of the table printed after the solve.
- **Reference cases** (U-AMG):
  - Cube-tet, k=0: `./bin/fhhos4 -geo cube -mesh tetra -k 0 -n 32 -s fcguamg -no-cache` (318k face unknowns)
  - Cube-cart-aniso100: `./bin/fhhos4 -geo cube -mesh cart -mesher inhouse -k 0 -n 64 -aniso 100 -s fcguamg -no-cache`
  - High order: `./bin/fhhos4 -geo cube -mesh tetra -k 2 -n 16 -s fcguamg -hp-cs h -no-cache` (262k face unknowns).
    `-hp-cs h`: the h-coarsening of the k=2 blocks, the default until 2026-10-05, on which the k=2 measurements of
    this file were made. The default is now p_h (see Done).
  - The variants of Fig. 4.5 / Table 4.4: `-coarsening-prolong 6` (P_F, default), `4` (Q_F^smooth), `5` (P_F^(0)), `3` (Q_F).
- **Profiling**: no sampling profiler, so temporary wall-clock timers around the steps, per multigrid level:
  `scripts/perf/Prof.h` (`PROF_START`/`PROF_STOP`, printed at exit). Include it in the files to instrument, never
  commit the instrumentation.
- **Validation of a change**:
  1. `ctest` (128 tests).
  2. `scripts/perf/amg_papers.sh <output dir> [binary]`: the AMG configurations of the papers, at reduced sizes (the
     paper sizes don't fit in 13 GB), 30 runs in ~6 minutes. Then `scripts/perf/fingerprint.sh <dir before> <dir after>`
     compares their iteration tables (iteration, inner iterations, residual): they must be identical.
     The runs: Cube-cart n=32, Cube-tet n=16, Complex-tet (`-geo platewith4holes`) n=4, Heterog1e8 n=256,
     Cube-cart-aniso100 n=32, Cube-tet-aniso20 n=16 and Cube-tet n=32 (Fig. 4.4), each with `fcguamg` and
     `fcgaggregamg`; Cube-tet n=16 with `-cs mpa|dpa`; Cube-tet n=16 and Cube-cart-aniso100 n=32 with
     `-coarsening-prolong 3|4|5|6`; the biharmonic Table 3 runs (n=8,16, `-bihar-prec s|no`) and the commented-out
     k=3 runs (`-opt2 0|2`) of `2023_biharmonic_problem.md`. After a change that reorders operations, the residuals
     may differ in their last digits: the iteration counts must be identical.
  3. When replacing a computation: temporary bit-for-bit checks of each new matrix against the old computation,
     with `scripts/perf/BitCheck.h` (same structure, `memcmp` of the values), on the runs of `amg_papers.sh`
     (`grep BITCHECK` in the logs).
  4. For a change of the solve phase (smoothers, cycle): `scripts/perf/solutions.sh <output dir> <binary before>
     [binary after]` compares the solution vectors (exported with 17 digits: exactly) and the iteration tables of 32
     cases covering the code paths of the smoothers, sequential and parallel, in ~10 minutes. A bit-identical change
     gives `OK` everywhere, except the known differences listed in the script. A change that reorders operations gives
     `DIFF` and `RESIDUALS-DIFF` (expected), but never `ITERATIONS-DIFF`.
- **Setup only (U-AMG)**: `scripts/perf/UAMGDump.h` dumps the inputs of `UncondensedAMG::Setup()` (the condensed
  operator, A_TT, A_TF, A_FF and the parameters of the multigrid; temporary instrumentation, see the header), and
  `scripts/perf/uamg_harness.cpp` runs the setup alone on them (`build_uamg_harness.sh`: ~30 s to compile with the
  flags of the build, against a ~2.5-minute rebuild of `Program.cpp`, and no assembly: 178 s sequential at k=2 n=16).
  Its setup times are within 5-15% of `timing.sh`'s. It can save the operator and the prolongation of each level,
  and compare two such runs (same structure? largest differences): the quick check of a change of the setup, before
  the validation above. Compile it against a `git archive` of HEAD for the "before" times.
- **Kernel micro-benchmarks**: the matrices of a run can be exported with `-export mg -o <dir>` (`levelN_A.dat`,
  `levelN_P.dat`, `levelN_R.dat`, Matrix Market, 17 digits) and loaded with `Eigen::loadMarket` in a standalone
  program compiled with the flags of the build (`grep FLAGS build/build.ninja`): iterating on a kernel then takes
  seconds instead of a 2-minute rebuild of `Program.cpp`. Check its outputs with `memcmp` against the old kernel.

## Exactness facts (to keep the results bit-identical)

- **Eigen's sparse product** `A * B` (row-major): each coefficient is summed in the storage order of the row of `A`,
  the first term is assigned (not added to 0), the explicit zeros are kept, and the rows are sorted.
  `SparseMatrixOps::Multiply` does the same, whatever the number of threads. Eigen evaluates the products from left
  to right: `P.transpose() * S * P` is `Multiply(Multiply(Pt, S), P)`.
- **`S.selfadjointView<Eigen::Lower>()`** only reads the lower triangle: `SparseMatrixOps::FullFromLower(S)`. It
  matters: the Galerkin operators are not bitwise symmetric.
- **`NonZeroCoefficients::Add`** drops the coefficients of absolute value <= `NonZeroCoefficients::ZeroThreshold`
  (1e-15). `setFromTriplets` sorts the rows and sums the duplicates.
- These facts matter when a change must be bit-identical (e.g. a refactoring, to validate it with `solutions.sh`).
- **Eigen's sparse matrix-vector product** (row-major, Eigen 3.5 of the conda env): each row is summed with two
  accumulators, for the even and odd positions of the row, added at the end (`res += alpha * (even + odd)`). Its
  sparse triangular solves (`triangularView<Lower|Upper>().solve()`) subtract the products one by one, in the order of
  the row, then divide by the diagonal. The build targets `-march=nocona` (no FMA): `a*b + c` is a multiplication then
  an addition, so hand-written loops reproduce these exactly (`GaussSeidel`'s kernels do).
- **Products by zero**: `s - a*0 = s` (exactly, up to the sign of a zero `s`): the products by an initial guess x = 0
  can be skipped (`BlockGaussSeidel` does). The flag `xEquals0`, passed down the cycle, says when x = 0: in the first
  sweep of the pre-smoother of every coarse level (the coarse correction starts from 0, `Multigrid.h`), and of the
  fine level when the multigrid is the preconditioner of FCG (`fcguamg` etc.: each application solves A e = r from
  e = 0, `Preconditioner.h`). Not on the fine level of a standalone multigrid (`-s uamg`, `mg`...) after its first
  iteration, nor in the later sweeps of a pre-smoother with several iterations, the second visit of a coarse level
  in a W-cycle, or the post-smoother.
- **The hybrid smoothers** (`hbgs`, `hrbgs`: Gauss-Seidel inside the chunk of rows of each thread, Jacobi between
  the chunks) give results that depend on the number of threads, and vary from one parallel run to the next: a
  thread reads x while another one writes it.

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

- **Smoothers** (2026-10-03): the default smoothers of U-AMG (`bgs`/`rbgs`), and of the other multigrids. Bit-identical
  (`solutions.sh`: all cases `OK`). The positions of the diagonal blocks in the rows are computed once
  (`DiagonalBlockPositions`); the matrix must be compressed, with its rows sorted by column index (checked; Eigen's
  matrices always are: `setFromTriplets`, `BuildByRows`).
  - `GaussSeidel` (block size 1, i.e. k=0): each sweep reads the rows once instead of twice (Eigen's `b - U*x` then
    triangular solve), the diagonal positions are computed in the setup, the residual/`Ax` loops are parallel (rows
    independent, like Eigen's product). Fine-level kernels (micro-benchmark, k=0 n=32): forward sweep from 0 +
    residual 6.2 -> 4.7 ms, backward sweep + `Ax` 9.5 -> 6.3 ms.
  - `BlockGaussSeidel` (k>=1): the sum of each row in a register (it was accumulated in an Eigen vector: a
    store/load chain per coefficient), no integer division per coefficient (`j / blockSize`: the diagonal block is
    located in the setup), in backward sweeps the rows of a block row processed together (independent chains, and
    their memory accesses start together; in forward sweeps, row by row is faster: one sequential stream), and the
    products by the zero initial guess skipped (when `xEquals0`, see the exactness facts). Fine level, k=2 n=16:
    micro-benchmark, forward sweep 28 -> 11 ms, backward sweep 37 -> 18 ms; in the solve (35 iterations, sequential),
    pre-smoother (forward sweep from 0 + residual) 1.49 -> 0.75 s, post-smoother (backward sweep + `Ax`) 1.49 -> 0.80 s.
    The hybrid smoother (`hbgs,hrbgs`, k=2 n=16, 16 threads): solve 3.2-3.6 -> 2.6 s.
  - `BlockSOR` (`bsor` etc. with omega != 1, and the standalone solvers `bgs`/`rbgs`): the same sums in registers,
    without division.

  Solve time (best of 2), before -> after (setup unchanged):

  | Case | Sequential | Parallel (16 threads) |
  |---|---|---|
  | Cube-tet k=0 n=32, P_F | 0.95 -> 0.81 s | 0.95 -> 0.81 s |
  | Cube-tet k=0 n=32, Q_F^smooth / P_F^(0) / Q_F | 0.94 / 1.04 / 1.04 -> 0.76 / 0.85 / 0.88 s | 1.00 / 1.01 / 1.08 -> 0.89 / 0.89 / 0.90 s |
  | Cube-cart-aniso100 n=64 | 0.59 -> 0.47 s | 0.61 -> 0.52 s |
  | Cube-tet k=2 n=16 | 6.24 -> 3.71 s | 6.95 -> 3.56 s |

  Tried, no gain: a hand-written solve of the small diagonal blocks (bit-identical to Eigen's LU solve): same time,
  the block solve is latency-bound (chain of divisions), not overhead-bound; software prefetch of the next block row
  in the backward sweep: up to -10% on the row-by-row version, still slower than processing the rows together; glibc's
  malloc thresholds (`MALLOC_MMAP_THRESHOLD_`, the cycle allocates its vectors): no effect; in the scalar sweeps,
  multiplying by the precomputed inverse of the diagonal instead of dividing by it (allowed, not bit-identical):
  0.92x-1.12x, noise (micro-benchmark, k=0 n=32, levels 0 and 1).

- **Hybrid meshes of U-AMG** (2026-10-04): bit-identical (temporary checks against the old computations on the runs
  of `amg_papers.sh`, plus k=1 and k=2 in 2D and 3D and `-prolong 5|6`: 13.9M couplings and 258 `Theta()`, no
  mismatch; `fingerprint.sh` identical; `ctest` 85/85).
  - The auxiliary coarse mesh of each pass, and the coarse mesh of the non-chained prolongations (`-prolong` != 2), only
    build the faces of their elements (`HybridAlgebraicMesh::BuildElementFaces()`): the prolongations only use those
    (`Theta()`) or the matrices (`-prolong 7`). No faces' elements, column-major copy of `A_T_F`, neighbours and
    couplings.
  - The couplings of `Build()` come from the traces of the blocks `A_T_F(element, face)`, summed once per element and
    face during the pass over the rows of `A_T_F`, in the order of `DenseMatrix::trace()`. They were computed from two
    dense blocks extracted with `block()` (heap allocations, scan of the rows) per neighbour and direction.

  Timers around the two steps (all passes and levels) and setup (`timing.sh`, best of 2), before -> after:

  | Case | `mesh.Build()`, seq. / 16 threads | Aux. coarse mesh, seq. / 16 threads | Setup, seq. / 16 threads |
  |---|---|---|---|
  | Cube-tet k=0 n=32, P_F | 0.59 -> 0.43 / 0.14 -> 0.14 s | 0.25 -> 0.05 / 0.11 -> 0.06 s | 4.29 -> 3.38 / 1.96 -> 1.74 s |
  | Cube-cart-aniso100 n=64 | 0.50 -> 0.46 / 0.20 -> 0.21 s | 0.30 -> 0.09 / 0.16 -> 0.11 s | 4.67 -> 4.32 / 2.73 -> 2.63 s |
  | Cube-tet k=2 n=16 | 0.33 -> 0.20 / 0.15 -> 0.12 s | 0.15 -> 0.03 / 0.11 -> 0.01 s | 8.20 -> 7.84 / 3.23 -> 3.16 s |

  Over the 4 coarsening prolongations at k=0, sequential: `mesh.Build()` 0.55 -> 0.45 s, auxiliary mesh 0.23 ->
  0.05 s, i.e. setup -0.2 to -0.36 s (-6 to -10%); on 16 threads, -0.01 to -0.09 s. The sequential setup before
  with P_F (4.29 s) is a slow outlier (its solve, unchanged, too: 0.89 vs 0.74 s): the timers give -0.36 s. The
  couplings were a smaller part of `Build()` than the profile suggested (0.39 s for its whole neighbours loop): the
  rest of that loop (vectors growing per element, sort, a mutex per element, scattered accesses to the neighbours)
  remains, see idea 2.

- **Local Galerkin products in U-AMG** (2026-10-04): not bit-identical (sums in another order), iteration counts
  unchanged. `src/Solver/Multigrid/UncondensedAMG/LocalMatrices.h`. A_TT is block-diagonal, so the condensed operator
  is a sum of local matrices on the faces of the cells, `S = sum_T S_T`, and P_F sends the faces of a cell only to the
  coarse faces of its aggregate (kept faces: injection; removed faces: smoothed reconstruction in the aggregate). So
  `P^T S P = sum_T P_T^T S_T P_T` (P_T: the rows of P on the faces of T), a sum of dense local matrices S_K on the
  aggregates. The coarsening passes go on with them, from one pass and one level to the next (`LocalOperator`, shared
  by the levels): the operator of a level is only assembled at its end, for the smoothers and the coarse solver. The
  block Jacobi of P_F uses the rows of the removed faces assembled from the local matrices, and only factorizes their
  diagonal blocks (`BlockJacobi::SetupForIterationMatrix()`, also used with the global products: bit-identical there).
  - The fine operator is decomposed once: `S_T(f,g) = S(f,g) / c`, where c is the number of cells containing both f
    and g. Exact for c <= 2: the reassembled S is bit-identical to `FullFromLower(S)` (Cube-tet k=0 and 2, Cube-cart,
    Square-tri k=3, polygonal meshes). Off the diagonal, c = 1 on the simplicial and Cartesian meshes; on the
    polygonal meshes built with `-polymesh-fcs n`, ~2400 pairs of faces are shared by 2 cells at n=64: no fallback
    needed. On the coarse levels, two coarse faces share at most one aggregate (one coarse face per pair of
    neighbouring aggregates), but the local products don't rely on it.
  - Used with `-face-prolong 1|2` (interface collapsing) and `-coarsening-prolong 3|4|5|6`; the other options keep the
    global products (with `-face-prolong 3`, the faces are aggregated independently of the cells: P leaves the
    aggregate).
  - The coarse operators are now exactly symmetric (the global products were not, at the rounding level).
  - Checks: on every pass of Cube-tet k=0 n=32 and k=2 n=16 and of 20 other configurations (Cube-cart, anisotropic,
    Complex-tet, Heterog, Square-cart, biharmonic, k=1-3, `stri`, `-coarsening-prolong 3|4|5`, `-cs dpa|mn`,
    `-face-prolong 2`, `-fcs i`, `-prolong 5`), the local product had the structure of the global one, and
    `|diff| / sqrt(|d_i d_j|) <= 1.2e-15`; per coefficient, relative differences up to 4e-9 on cancelled coefficients,
    as large as the global product's own asymmetry. Whole setup (`uamg_harness`): the operator and the prolongation of
    every level have the same structure as before (same aggregates), `|diff| / sqrt(|d_i d_j|) <= 8e-15`. `ctest`
    92/92; `amg_papers.sh`: identical iteration tables (even the 3-digit residuals); `solutions.sh`: `DIFF`,
    `RESIDUALS-DIFF` in 2 cases, and one `ITERATIONS-DIFF`: the hybrid smoother in parallel, which gives 37 or 38
    iterations from one run to the next, before too. Identical iteration counts also with `-fcs i`, `-cs mn`,
    `-face-prolong 2`, `-prolong 5|6`, `-hp-cs p_h|h_p`, Cube-cart k=1, Square-tri k=1, `stri` k=3, Heterog k=1,
    `poly` (`-polymesh-fcs c|n`) and the in-house `quad` mesh (`-face-prolong 3`, which keeps the global products, runs
    for more than 10 minutes on Cube-tet k=0 n=16, before as after: not compared).
  - Memory: the local matrices of a pass take about as much as S (k=2 n=16: 105 MB for the cells of the first pass).
    Two buffers, reused from one pass to the next (see the first access to fresh memory, in Machine), freed at the end
    of the setup.

  The replaced operations took (all passes, sequential / 16 threads) 4.8 / 1.55 s at k=2 n=16 (`Pt*S` 1.99 / 0.42 s,
  `(Pt*S)*P` 1.52 / 0.33 s, `FullFromLower(S)` 1.27 / 0.75 s, transposes 0.03 / 0.04 s) and 0.89 / 0.35 s at k=0 n=32.
  Their replacements take 0.56 / 0.23 s at k=2 (products 0.29 / 0.12 s, decomposition of the fine operator 0.12 /
  0.04 s, rows of the removed faces 0.07 / 0.03 s, assembly of the operators of the levels 0.06 / 0.03 s, block
  Jacobi factorizations 0.02 / 0.01 s) and 0.59 / 0.20 s at k=0 (products 0.28 s, rows of the removed faces 0.15 s,
  decomposition 0.11 s, assembly 0.04 s sequential).

  Setup (`timing.sh`, best of 2), before -> after; the solve is unchanged:

  | Case | Sequential | Parallel (16 threads) |
  |---|---|---|
  | Cube-tet k=0 n=32, P_F | 3.49 -> 3.20 s | 1.96 -> 1.71 s |
  | Cube-tet k=0 n=32, Q_F^smooth / P_F^(0) / Q_F | 3.31 / 3.19 / 3.05 -> 2.95 / 2.82 / 2.62 s | 1.91 / 1.87 / 1.76 -> 1.64 / 1.62 / 1.50 s |
  | Cube-cart-aniso100 n=64 | 4.45 -> 4.20 s | 2.63 -> 2.49 s |
  | Cube-tet k=2 n=16 | 8.07 -> 3.68 s | 3.39 -> 1.78 s |

  Cube-tet k=0 n=48 (16 threads, 2 runs each): setup 7.1-7.7 -> 6.2 s, solve unchanged (3.0 s), peak memory of the
  whole run 4.32 -> 4.04 GB. At k=0, `P^T S P` was ~18% of the setup on 16 threads (~45% at k=2): the hybrid meshes
  are now the largest part (idea 2).

  How it was developed: a prototype computing both products in each pass (`uamg_harness`), then the decomposition
  of the fine operator only, and the assembly only at the end of a level. A first prototype was 2.4x faster than the
  global products at k=2 but slower at k=0: a heap allocation per cell (1x1 to 6x6 blocks at k=0), a sort per row in
  the assembly, and the first access to ~400 MB of fresh memory per pass (page faults); reused buffers, rows written
  in place (their lengths counted first) and plain loops instead of Eigen's dynamic-size products fixed it.

- **U-AMG at k >= 1: p-levels, then h, by default** (2026-10-05): `-hp-cs p_h` instead of `h` when k >= 1
  (`ProgramArgumentsDefaults.cpp`). A change of the method, decided by the user: the h-coarsening of the degree-k
  blocks transfers their higher modes by plain aggregation (identity blocks of Q_T, Q_F; trace of the constant
  only), whereas p_h only h-coarsens at k=0, where the operators of the paper apply. The published U-AMG runs are
  at k=0: `amg_papers.sh` gives identical iteration tables, except the commented-out k=3 biharmonic runs (outer
  iterations 27 -> 27 and 35 -> 31, solve 6.6 -> 2.9 s and 3.6 -> 1.4 s with `-opt2 0` and `2`).
  32 cases (`fcguamg`, 16 threads, one run each, logs in `build/bench_hp/`), setup + solve:
  - Cartesian, triangular, polygonal, heterogeneous (2D/3D, k=1..3): 1.5-4x fewer iterations, 1.2-4x faster.
  - Tetrahedral (Cube-tet, Complex-tet, Cube-tet-aniso20): same iterations (the k=0 U-AMG sets them), 1.2-3x
    faster: operator complexity 1.5-1.7 -> 1.04-1.2. Cube-tet k=2 n=16: setup 1.8 -> 0.3 s, solve 3.4 -> 2.1 s.
  - Anisotropic Cartesian (aniso100, 2D/3D): 0-5 more iterations, as fast or faster (cheaper cycle).
  - Peak memory -10 to -25% (Cube-tet k=1 n=32: 4.27 -> 3.28 GB). Same growth of the iterations with n.

  The high-order reference case of `timing.sh` keeps `-hp-cs h`, so that its series stays comparable: the k=2
  measurements of this file (profile, ideas) concern that path, no longer the default.

- **Library `fhhos4_AMG`** (2026-10-07, `library/README.md`): no change of the results of the program (`amg_papers.sh`:
  identical iteration tables, 30 runs; `ctest`). Changes of its paths: `Utils::FatalError()` throws instead of
  `exit()` (no effect when nothing is thrown: not measured); the block `A_F_F` is optional and no longer extracted on the
  p-levels when the options don't use it (default: the algorithm only uses `A_T_T` and `A_T_F`), which a client then
  doesn't need to assemble. In HArDCore3D (`hho-diffusion -c 1 2`, sequential, tolerance 1e-12, solve wall time
  including the client's assembly of `A_TT` and `A_TF`): 32^3 k=1 2.5 s (FCG + K-cycle, 19 iterations) vs 4.6 s for its
  Jacobi-BiCGSTAB (260 iterations), k=2 8.3 vs 17.3 s, 48^3 k=1 9.5 vs 25.5 s; BiCGSTAB + V-cycle 4.0 / 13.3 / 18.4 s.
  Beware of HArDCore3D's test case 1 (`sin(pi x) sin(pi y) sin(pi z)`, constant diffusion): the first eigenfunction of
  the Laplacian, so a Krylov method converges in a few iterations, whatever the preconditioner.

- **Library: hybrid Gauss-Seidel by default on several threads** (2026-10-08, idea 6). The smoothers of
  `fhhos4::Solver` default to empty: `Setup()` then chooses `hbgs`/`hrbgs` (Gauss-Seidel in the rows of each thread,
  Jacobi between the threads) when `Krylov` is `fcg` and several threads run, `bgs`/`rbgs` otherwise. The program is
  unchanged (`-smoothers bgs,rbgs`: the papers), and so are its library runs (`-s libuamg`, which pass the program's
  smoothers). Done in two steps: at k >= 1 first, with the existing hybrid smoother; then at k = 0 too, once
  `hbgs`/`hrbgs` with blocks of size 1 ran the scalar kernels of `GaussSeidel` (its hybrid mode, `HybridSweep()`:
  the sequential kernels on the chunk of each thread; the residual computed apart, since the fused residual of the
  sequential sweeps doesn't hold between chunks; from x = 0, the products by the values not computed yet are skipped,
  as in the sequential sweeps; parallel above 20000 non-zeros). Before, they ran `BlockGaussSeidel`'s code with 1x1
  blocks. The default sequential sweeps are unchanged (identical residual tables, `-threads 1`).

  Solve wall time, `bgs,rbgs` -> `hbgs,hrbgs`, FCG + K-cycle (setup unchanged, 1-2 runs each). HArDCore3D:
  `hho-diffusion -c 1 2 --solver_type fhhos4`, tolerance 1e-12, through a temporary override of the smoothers in the
  library; fhhos4: `-s fcguamg`, tolerance 1e-8. First step (block code at k = 0, i.e. on the k=0 levels too):

  | Case | 1 thread (`bgs`) | 8 threads | 16 threads |
  |---|---|---|---|
  | HArDCore3D 32^3 k=1 (19 its) | 1.08-1.10 s | 1.06-1.17 -> 0.72-0.74 s | 1.33-1.43 -> 0.97-1.22 s (20 its) |
  | HArDCore3D 32^3 k=2 (26 its) | 4.61 s | 4.00-4.11 -> 2.44-2.75 s | 4.78-4.82 -> 2.95-3.66 s |
  | HArDCore3D 48^3 k=1 (20 its) | 4.94 s | 4.19-4.59 -> 2.58 s | 4.17-4.22 -> 2.83-2.90 s |
  | HArDCore3D voro-16 k=1 (16 its) | 0.73 s | 0.71-0.72 -> 0.48-0.49 s | |
  | HArDCore3D voro-16 k=2 (20 -> 21 its) | 2.89 s | 2.71-2.76 -> 1.81-1.82 s | |
  | fhhos4 Cube-tet k=1 n=16 (33 -> 38 its) | | 0.59-0.60 -> 0.46-0.47 s | 1.00-1.31 -> 0.65-0.70 s (37-39 its) |
  | fhhos4 Cube-tet k=2 n=16 (33 -> 37-38 its) | | 1.77-1.85 -> 1.20-1.23 s | 2.01-2.10 -> 2.39-2.53 s (39 its) |

  On 2 threads, HArDCore3D 32^3 k=2: 3.86 -> 3.15 s (26 its). The tetrahedral meshes of GMSH lose more iterations
  (their face numbering has little locality: more couplings between the chunks of the threads). L2 errors of
  HArDCore3D identical to the references (`bench_lib/quick_check.sh`, and its 4 runs on 8 threads). Not chosen:
  - k = 0: the hybrid smoother runs the block code with 1x1 blocks, slower than the scalar `GaussSeidel`. HArDCore3D
    48^3 k=0, 8 threads: 0.77 -> 0.86-0.94 s (32 -> 33 its); fhhos4 Cube-tet k=0 n=32: 16 threads 0.85-1.04 ->
    0.97-1.30 s (28 -> 31 its), 1 thread 0.75 -> 1.64 s.
  - 1 thread: the same sweeps in slower code (the products by the zero initial guess are not skipped, 1x1 blocks on
    the k=0 levels): HArDCore3D 32^3 k=1 1.08 -> 1.45 s, voro-16 k=1 0.73 -> 0.93 s, 48^3 k=1 4.94 -> 5.50 s.
  - BiCGSTAB + V-cycle: the hybrid smoother varies from one application to the next (timing of the threads), which
    BiCGSTAB doesn't expect. 8 threads, 32^3 k=1: 26 -> 37-38 its (2.37 -> 1.78-1.83 s); k=2: 31 -> 49 its (8.03 ->
    6.51 s). Faster, but not as a default.

  Tried, no gain: hybrid on the levels with blocks, scalar sequential `GaussSeidel` on the k=0 levels of the p-then-h
  hierarchy (8 threads: 32^3 k=1 0.68-0.70 s, k=2 2.28-2.33 s, 48^3 k=1 2.54-2.86 s, voro-16 k=1 0.53-0.57 s: within
  the noise of the plain hybrid smoother).

  Second step, with the scalar hybrid kernels at k = 0, 8 threads (the new default, see the next entry):

  | Case | `bgs,rbgs` | `hbgs,hrbgs` |
  |---|---|---|
  | HArDCore3D 32^3 k=0 | 0.20 s (29 its) | 0.155-0.170 s (30-31 its) |
  | HArDCore3D 48^3 k=0 | 0.75-0.78 s (32 its) | 0.63-0.64 s (33 its) (block code: 0.86-0.94 s) |
  | HArDCore3D voro-16 k=0 | 0.090-0.095 s (21 its) | 0.074-0.077 s (21 its) |
  | fhhos4 Cube-tet k=0 n=32 | 0.59 s (28 its) | 0.59-0.60 s (33-34 its) |
  | fhhos4 Cube-cart-aniso100 n=64 | 0.42-0.43 s (9 its) | 0.44-0.46 s (11 its) |
  | fhhos4 Square-tri k=0 n=512 (mesh from the cache) | 2.03 s (26 its) | 1.71-1.75 s (29 its) |
  | HArDCore3D 32^3 k=1 | 1.06-1.10 s (19 its) | 0.61-0.65 s (19 its) (block code on the k=0 levels: 0.72-0.74 s) |
  | HArDCore3D 32^3 k=2 | 4.00-4.11 s (26 its) | 2.29-2.35 s (26 its) (2.44-2.75 s) |
  | HArDCore3D 48^3 k=1 | 4.19-4.59 s (20 its) | 2.29-2.32 s (20 its) (2.58 s) |
  | HArDCore3D voro-16 k=1 / k=2 | 0.71 / 2.71-2.76 s (16 / 20 its) | 0.43-0.47 / 1.73 s (16 / 21 its) (0.48 / 1.81 s) |

  Hence the default at k = 0 too: faster on HArDCore3D's meshes and on the 2D triangles, as fast on fhhos4's
  tetrahedral and anisotropic meshes. Still not chosen on 1 thread (the hybrid sweeps then skip fewer products and
  don't give the residual). L2 errors of HArDCore3D identical to the references and to its BiCGSTAB (k = 0..2, 1 and 8
  threads).

  Then whatever the Krylov method (the user's request: one default): with BiCGSTAB, which expects the same
  preconditioner at each iteration, the hybrid smoother takes more iterations, but it is faster with every method
  (HArDCore3D, 8 threads, `bgs,rbgs` -> `hbgs,hrbgs`): BiCGSTAB + V-cycle, 48^3 k=0 1.13 s (29 its) -> 1.04-1.11 s
  (39-40 its), 32^3 k=1 2.55 s (26 its) -> 1.91 s (38 its); the multigrid alone (`Krylov = "none"`, K-cycle), 48^3
  k=0 3.34 s (119 its) -> 2.18-2.23 s (124 its), 32^3 k=1 1.86 s (31 its) -> 1.02 s (32 its).

- **The solver on one thread per physical core by default** (2026-10-08, the user's decision): the solvers run on at
  most one thread per physical core (no hyper-threading), the rest of the program (mesh, assembly, post-processing)
  on all the logical cores, the OpenMP default. In the program, the solver section of each program (creation, setup,
  solve) is a `Parallelism::SolverThreads` scope; in the library, `Threads = 0` does the same during each call.
  `Parallelism::WithoutHyperThreads()`: no cap if `OMP_NUM_THREADS` is set, nor with OpenMP's binding
  (`OMP_PROC_BIND`, `OMP_PLACES`: libgomp then pins the main thread to its place at startup, so its affinity no longer
  gives the cores of the process: checked), nor if the cores are unknown; an explicit `-threads N` applies to the
  solver too. `Parallelism::PhysicalCores()`: the number of distinct `thread_siblings_list` among the CPUs of the
  affinity mask (`/sys`, Linux; elsewhere 0, no cap), so it follows `taskset`, Slurm and MPI bindings (checked: 3 cores
  for `taskset -c 0-5`, 4 for `taskset -c 0,2,4,6`). No portable API gives the physical cores (hwloc would, as a
  dependency). Results unchanged (the default solvers don't depend on the number of threads; `ctest` 130/130).
  Measured under WSL2 on this machine (on bare-metal Linux, the hyper-threading penalty of memory-bound code is often
  smaller), 16 -> 8 threads:
  - Solver, `bgs` (the program's default smoother): HArDCore3D 32^3 k=2, setup 1.8-2.1 -> 0.48-0.77 s, solve 4.8 ->
    4.0-4.1 s; 32^3 k=1, setup 0.51-0.86 -> 0.41-0.45 s, solve 1.33-1.43 -> 1.06-1.17 s (slower than on 1 thread,
    1.08 s, on 16 threads); fhhos4 Cube-tet k=0 n=32, solve 0.85-1.04 -> 0.59 s; Cube-tet k=1 n=16, solve 1.00-1.31 ->
    0.59-0.60 s; k=2 n=16, 2.01-2.10 -> 1.77-1.85 s.
  - Solver, `hbgs` (the library's default): HArDCore3D 32^3 k=1, setup 0.57-0.61 -> 0.43-0.48 s, solve 0.91 -> 0.68 s
    (20 -> 19 its); 32^3 k=2, setup 0.62-0.77 -> 0.55-0.61 s, solve 3.18-3.41 -> 2.54-3.14 s; 48^3 k=0, solve 1.40-1.54
    -> 0.72-0.91 s (setup 1.5-2.0 s both); 48^3 k=1, setup 1.96-2.83 -> 1.63 s, solve 2.99-3.18 -> 2.55 s.
  - The assembly prefers all the logical cores where it is compute-bound, hence the split: Cube-cart-aniso100 n=64
    (the quadrature of the right-hand side, see Current profile), 8 -> 16 threads, 26.9-29.2 -> 21.8-24.1 s in one
    series, 21.8 -> 13.9 s in another (the machine's speed varies: CPU time per run +50% from one series to the next);
    Cube-tet k=2 n=16, 44.7 s on 8 threads, 38.2 s on 16, 36.4 s by default (solve 1.92 s on 8, 2.31 s on 16, 1.85 s by
    default).
  `OMP_WAIT_POLICY=passive` (idea 1) was re-measured with these runs: slower.

- **Library: row-major matrices read in place** (2026-10-09). `Setup()` copied every input matrix into fhhos4's format
  (row-major, `int` indices, compressed), and kept the copy of `A` (the levels point to it). The matrices already in
  this format are now read in place, and the solver keeps a reference to the caller's `A` (it must outlive the solver;
  a temporary `A` doesn't compile); the others (column-major, `long`, not compressed) are still converted. Results
  unchanged (`ctest` 131/131, HArDCore3D quick check: same L2 errors and iterations). The program's `-s fcglibuamg`
  (fhhos4's matrices), Cube-cart k=1 n=32, 3 runs each: setup 0.348-0.365 -> 0.267-0.284 s, peak memory 941-943 ->
  875-879 MB, as `-s fcguamg` (0.295 s, 872 MB). Unchanged for HArDCore3D, whose matrices are column-major.

## Current profile (U-AMG; parallel: 16 threads, the default before 2026-10-08)

| Case | Setup, sequential / parallel | Solve, sequential / parallel |
|---|---|---|
| Cube-tet k=0 n=32 | 3.2 / 1.7 s | 0.7-0.9 / 0.7-0.9 s |
| Cube-cart-aniso100 n=64 | 4.2 / 2.5 s | 0.5 / 0.5 s |
| Cube-tet k=2 n=16 | 3.7 / 1.8 s | 3.7 / 3.5 s |

**The assembly (printed as "Assembly time", not in the Setup row) dominates the runs** (2026-10-04, `timing.sh`
logs, sequential / 16 threads): Cube-tet k=2 n=16 178 / 38 s (setup 3.7 / 1.8 s); Cube-cart-aniso100 n=64 102 / 16 s
(setup 4.2 / 2.5 s); Cube-tet k=0 n=32 5.1 / 1.4 s. Breakdown (timestamps of the assembly log, n halved):
- The static condensation (`A_T_ndF^T * Solve_A_T_T(A_T_ndF)`, global product) takes 0.08 s: negligible.
- k=2 (Cube-tet n=8, 27.5 s): computing the local matrices takes 26.9 s (~9 ms per element). Not profiled yet.
- Cartesian meshes (Cube-cart n=32, 13.2 s, with or without `-aniso`): the right-hand side takes 12.5 s (0.1 s on a
  tetrahedral mesh). On a Cartesian cell, `ReferenceCartesianShape::Integral(RefFunction)` integrates a
  non-polynomial function with `GaussLegendre::MAX_POINTS` = 20 points per direction (8000 in 3D), and
  `PhysicalShape::InnerProductWithBasis` evaluates f again for each basis function; tetrahedra use Keast's default
  rule (degree 8). Fewer points change the right-hand side at the level of the quadrature error (opt-in, or a
  decision); evaluating f once per point for all the basis functions is bit-identical (10x fewer evaluations at k=2).
- Basis ids (2026-10-07: the matrices stored by the reference shapes are keyed by `FunctionalBasis::Id`, an atomic
  counter, instead of the address of the basis): assembly, setup and solve times unchanged on Cube-tet k=0 n=32,
  k=2 n=16, Cube-cart k=2 n=24 and the cases of `timing.sh` (Cube-tet k=2 n=16, 16 threads, interleaved runs:
  44.8 / 47.7 s before, 44.8 / 47.4 s after). Square-tri got faster: k=3 n=128 sequential 49-60 -> 39-43 s,
  k=1 n=256 10.2 -> 9.7 s. Not explained (no work removed; the bases are 8 bytes larger: memory layout?).

The breakdowns below, except the first, were measured before the optimizations of the hybrid meshes and of the
Galerkin products (see Done).

- **Setup after the local Galerkin products** (2026-10-04, `uamg_harness` with timers, sequential / 16 threads).
  k=2 n=16 (3.8 / 1.8 s): the coarsening passes take 2.7 / 1.3 s: building P_F 1.05 / 0.46 s (`Theta()` 0.53 /
  0.29 s, `BlockJacobi::IterationMatrix()` 0.21 / 0.05 s, `J * Y` 0.11 / 0.05 s, row selections 0.09 / 0.05 s), the
  local Galerkin products and their companions 0.56 / 0.23 s (see Done), the hybrid meshes 0.28 / 0.22 s; the rest
  (chaining of P, `R = P^T`, smoothers, coarse solver) 1.1 / 0.5 s. k=0 n=32 (2.9 / 1.7 s): the passes take 2.7 /
  1.3 s: hybrid meshes 1.17 / 0.67 s, building P_F 0.67 / 0.25 s, local Galerkin products and companions 0.59 / 0.20 s.
  See the ideas 2 to 5.

- **The solve barely benefits from the threads** (except with the library's default hybrid smoother, see Done).
  The sweeps of the default smoother (lexicographic Gauss-Seidel, as in the paper) are sequential; only the residual/`Ax` products and the intergrid transfers are parallel. Those are
  Eigen's: its row-major sparse matrix x vector product is already parallel (OpenMP, from 20000 non-zeros, dynamic
  schedule), so rewriting it adds no parallelism. And it is memory-bound: one thread already reads ~17.5 GB/s of the
  machine's ~35 GB/s, so the threads give at most ~2x on it (k=2, fine A: 7.4 -> 4.1 ms). Eigen's sequential
  operations are the sparse x sparse products and the transposes (in the setup: replaced by `SparseMatrixOps`), the
  products by a transposed matrix or a `selfadjointView` (scatters; not in the default solve path), and the
  triangular solves (sequential by nature).
- **Solve, k=0, sequential (0.8-0.9 s)**: the fine level is 2/3 of it: smoothing 0.42 s (post-smoother, backward
  sweep + `Ax`: 0.24 s; pre-smoother, forward sweep from 0 + residual: 0.18 s), prolongation 0.11 s, restriction
  0.10 s; level 1: 0.15 s. The fine kernels are limited by their scattered accesses to x: the numbering of the fine
  faces (GMSH) has little locality (mean |i-j| ~38000 over the non-zeros, against ~1100 on level 1, numbered by the
  aggregation), the fine SpMV runs at ~9 GB/s against ~19 GB/s on level 1. The division by the diagonal is not
  what bounds the sweeps: multiplying by its precomputed inverse instead gives no measurable gain.
- **Solve, k=2, sequential (3.7 s)**: by level, L0 2.1 s, L1 1.0 s, L2 0.2 s. Smoothing 2.75 s, of which the
  residual/`Ax` products after the block sweeps take ~0.6 s; restriction 0.32 s and prolongation 0.30 s (P has 37
  non-zeros per row at k=2: about as many as A); coarse solve 0.01 s.
- **Setup, k=0, 16 threads (2.1 s)**: hybrid meshes 0.9 s (build 0.16 s + auxiliary coarse mesh 0.12 s + coarsening
  0.45 s, including the sequential pairwise aggregation 0.18 s and interface collapsing 0.19 s + freeing 0.17 s);
  sparse products 0.4 s; `FullFromLower(S)` 0.16 s; `BuildQ_F` 0.10 s; `Theta()` 0.08 s; other transposes 0.09 s.
- **Setup, k=2, 16 threads (3.7 s)**: sequential transposes 1.2 s (`FullFromLower(S)` 0.81 s, `R = P^T` 0.24 s,
  column-major copy of `A_T_F` in `HybridAlgebraicMesh::Build()` 0.20 s); products 1.2 s (including `P^T*S` 0.52 s
  and `(P^T*S)*P` 0.43 s); `Theta()` 0.28 s; P chaining 0.20 s; hybrid meshes 0.4 s.
- **The setup is a sequence of pairwise passes** (measured 2026-10-03, timers per pass). `-cs mpa` (default) repeats
  the pairwise aggregation until the face coarsening factor reaches 3.8: Cube-tet k=0 n=32, 13 passes over 4 levels
  (4 on level 0: the face count only drops by 1.45x per pass at first); Cube-cart-aniso100, 2 per level. Each pass
  rebuilds a hybrid mesh, builds P_F and computes `P^T S P`, and the non-zeros of S barely drop from one pass to the
  next (level 0, k=0: 2.18M, 2.18M, 1.91M, 1.60M; k=2: 10.6M, 10.3M, 8.9M, 7.4M), so every pass costs about as much as
  the first: the passes after the first take 46-48% of the setup at k=0, 27-38% on Cube-cart-aniso100 n=64, 53-57% at
  k=2. Sequential / 16 threads, over all passes: the pairwise algorithm itself (`PairwiseAggregation::Perform`) 0.14 /
  0.18 s at k=0 (4% / 9% of the setup), 0.01 s at k=2; the hybrid-mesh bookkeeping around it 1.66 / 0.91 s at k=0
  (43-44%), 0.62 / 0.51 s at k=2; the algebra (P_F, Galerkin products, chaining) 1.77 / 0.73 s at k=0, 7.5 / 2.8 s at
  k=2 (sequential: `P^T*S` 2.0 s, `(P^T*S)*P` 1.6 s, `FullFromLower(S)` 1.5 s; replaced since by the local Galerkin
  products, see Done).
- **Eigen's own kernels are at the hardware limit: neither rewriting them nor BLAS/LAPACK would help** (micro-benchmark,
  2026-10-03, exported matrices of Cube-tet k=0 n=32 and k=2 n=16, flags of the build). Sparse matrix-vector product
  (`A`, `P`, `R` of levels 0-1): hand-written CSR loops (1 or 4 accumulators) take the same time as Eigen's, within the
  noise. At k=2, Eigen's reaches 18 GB/s sequential and 32 GB/s on 16 threads, for a read bandwidth of 17.5 / 35 GB/s.
  At k=0, the fine `A` (12 GB/s) and `R` stay below it because of their scattered accesses to x (see above), which no
  SpMV implementation avoids. OpenBLAS (in the conda env; it selects AVX2 kernels at runtime): `ddot`/`daxpy` take the
  same time as Eigen (n ~300k, ~0.05-0.1 ms). On the 6x6 diagonal blocks (k=2), LAPACK `dgetrs` is 1.8x slower than
  Eigen's `PartialPivLU::solve` (126 vs 70 ns per block) and BLAS `dgemv` 2.1x slower than Eigen's product (43 vs
  20 ns): the call overhead dominates at this size. The rewritten smoothers are faster because they read the matrix
  fewer times, skip work and run in parallel (see Done), not because their primitives are faster than Eigen's. The
  remaining gains are in fewer bytes per non-zero (block CSR), locality, and fused passes.

## Next ideas, by expected gain

1. **`OMP_WAIT_POLICY=passive`** — measured, no code, results unchanged. k=2 (16 threads): solve 7.0-7.3 -> 6.2-6.3 s
   (-12%), setup unchanged or slightly better; k=0: no notable change. (`OMP_WAIT_POLICY=active` is much worse:
   solve 8.9 s.) libgomp reads it at startup, so it cannot be set from `main()`: set it in the environment, e.g.
   `conda env config vars set OMP_WAIT_POLICY=passive -n fhhos4`, and document it. Measured before the optimization of
   the smoothers, after which 16 threads are no longer slower than 1 at k=2 (3.6 vs 3.7 s): re-measure. Re-measured
   through the library (2026-10-08, HArDCore3D 32^3 k=1, solve): slower. 16 threads, `bgs` 1.33-1.43 -> 1.34-1.44 s,
   `hbgs` 0.97-1.22 -> 1.50-1.53 s; 8 threads, `bgs` 1.06-1.17 -> 1.20-1.22 s, `hbgs` 0.72-0.74 -> 1.08-1.11 s.
   Probably dropped: the library can't set it anyway (read at the start of the process).
2. **Flat data structures for the hybrid algebraic meshes** (`HybridAlgebraicMesh`) — estimated, medium-large
   effort, main lever at k=0 (the paper's case). Vectors of vectors, a mutex per element, a `std::map` of neighbours
   per aggregate, copied vectors: building and coarsening the meshes (`mesh.Build()` + `mesh.Coarsen()`) take 1.17 /
   0.67 s of the 2.9 / 1.7 s of setup at k=0, sequential / 16 threads (`uamg_harness`, after the local Galerkin
   products; 0.28 / 0.22 s of 3.8 / 1.8 s at k=2), plus freeing them. CSR-like arrays and prefix sums (same numbering
   of the aggregates and coarse faces, so bit-identical) might halve it. The pairwise aggregation itself (greedy, priority order) must stay sequential:
   parallelizing it changes the aggregates. Measured breakdown, over all passes, k=0 n=32, sequential / 16 threads
   (2026-10-03, before the optimization of 2026-10-04 in Done, which made the auxiliary meshes 0.26 -> 0.05 s and
   `Build()` 0.55 -> 0.45 s sequential): neighbours loop of `Build()` 0.39 / 0.09 s, of which the extraction of the
   couplings was ~0.1 s; auxiliary coarse meshes 0.26 / 0.18 s; aggregates' neighbours (`std::map`/`set`) 0.23 / 0.04 s;
   `BuildQ_T`/`BuildQ_F` (triplets, 1 non-zero per row for most rows) 0.21 / 0.12 s; removed faces 0.17 / 0.04 s;
   interface collapsing (sequential) 0.16 / 0.19 s; freeing 0.10 / 0.16 s; other `Build()` steps 0.13 / 0.10 s.
3. **`Theta()` without `SparseMatrix::block()`** — measured, small effort, bit-identical if the blocks and the solves
   are the same. `HybridAlgebraicMesh::Theta()` (the reconstruction of P_F, `-A_TcTc^-1 A_TcFc` per coarse cell)
   extracts its dense blocks with `block()`, which scans the row from its start for each face: 0.53 / 0.29 s of the
   3.8 / 1.8 s of setup at k=2 (sequential / 16 threads: it scales poorly), 0.15 / 0.07 s at k=0 (`uamg_harness`,
   2026-10-04). One pass over the rows of the cell, as in `BuildElementFaces()`, gives all the blocks; then a
   parallel fill without triplets. `CouplingValue()` (`-fcs i`) and `BuildHighOrderTraceOnRemovedFaces()` do the same.
4. **Local P_F** — estimated, medium effort, rounding-level differences. The rows of P_F on the removed faces only
   involve the faces of their aggregate: `P(f,:) = (1-w) Y(f,:) - w D_f^-1 sum_{g!=f} S(f,g) Y(g,:)` (block Jacobi,
   w = 2/3), with Y = the first row of Theta of the aggregate on the removed faces, the injection on the kept faces.
   Computed per aggregate from the local matrices (`LocalMatrices`), they would replace the assembly of the rows of
   the removed faces, `BlockJacobi::IterationMatrix()`, the product `J * Y` and the two row selections: 0.45 s
   sequential at k=2 (J 0.21 s, `J * Y` 0.11 s, rows of the removed faces 0.07 s, selections 0.06 s), 0.39 s at k=0
   (0.14, 0.07, 0.15, 0.03 s); ~0.15 / 0.12 s on 16 threads. They are also the rows the local Galerkin product reads
   in P. The assembled P is still needed (chaining, `A_T_Fc = Q_T^T A_T_F P`).
5. **Parallel transpose** — estimated, small effort, bit-identical (a transpose is exact). Replaces Eigen's sequential
   transposes in `UncondensedLevel::SetupRestriction` (`R = P^T`, 0.24 s on 16 threads at k=2 in the old profile),
   `HybridAlgebraicMesh::Build()` (column-major copy of `A_T_F`, 0.20 s), `SparseMatrixOps::Transpose` and
   `FullFromLower`. Since the local Galerkin products, the last two are no longer in the default path (only with the
   coarse `A_F_F` and the options that keep the global products). Re-measure: at most ~0.45 s of the 1.8 s of setup
   at k=2 on 16 threads, less at k=0.
6. **Parallel smoother** — the hybrid Gauss-Seidel (`hbgs,hrbgs`, scalar kernels at k = 0) is now the library's
   default on several threads (see Done: solve -15 to -45% on HArDCore3D's meshes); the program keeps `bgs,rbgs` (it
   changes the iteration counts and depends on the number of threads: opt-in only). What remains:
   - fhhos4's tetrahedral meshes (GMSH numbering, little locality): the hybrid smoother loses 5-6 iterations at
     k = 0 (28 -> 33-34) and gains nothing there. A multicolour Gauss-Seidel, or a renumbering of the chunks
     (see "Renumbering the fine unknowns"), might do better; gain unknown.
   - 1 thread: `BlockGaussSeidel::HybridSweep` could fall back to `Sweep` (bit-identical up to the sign of a zero,
     with the products by the zero initial guess skipped), and the scalar one could fuse the residual: `hbgs` would
     then cost nothing when a client runs it on 1 thread.
   - Determinism: the hybrid smoothers have a data race (see the exactness facts): a real Jacobi between the chunks
     (or a colouring of the chunks' boundaries) would make them deterministic, and usable with BiCGSTAB.

Possible now that reordering the arithmetic operations is allowed (rejected so far because not bit-identical):

- **Shorter dependency chains in the sweeps** — estimated, small effort. In `BlockRowRhs` (forward block sweeps,
  row by row) the sum of a row is one chain of ~40 dependent subtractions at k=2 (3-4 cycles each): 2 to 4 partial
  sums would bring the forward sweep (~11 ms on the fine level, k=2 n=16) closer to the SpMV (~7-8 ms).
- **Residual of the block pre-smoother from the sweep** — estimated, small effort. After a forward sweep,
  `r_I = U_I (x_old - x_new)` (in exact arithmetic, since `D_I x_I(new) = b_I - L_I x(new) - U_I x(old)`); from
  x = 0, `r_I = -U_I x(new)`: half of the matrix instead of the full `b - A*x` (~8 ms per pre-smoothing on the fine
  level at k=2, i.e. ~0.15 s of the 3.7 s solve). The scalar `GaussSeidel` already computes its residual this way.
- **Block CSR storage for k>=1** — estimated, large effort. One column index per block instead of one per
  coefficient (12 -> ~8 bytes per non-zero at k=2) and dense block-vector products (vectorizable): the solve kernels
  are memory-bound, so maybe -20-30% on the sweeps, residuals and intergrid transfers at k>=1.
- **`-DENABLE_NATIVE_ARCH=ON`** (AVX2, FMA) — to measure. FMA and wider vectors change the rounding (now allowed);
  the binary then only runs on CPUs like the build machine's.
- **Inverse diagonal blocks** in the block smoothers (`DiagBlockSolveMethod::Inverse`, an existing option of
  `BlockDiagonalSolver`, not the default): the solve of a 6x6 block becomes a matrix-vector product, without the
  chain of divisions of the LU solve. Measured (micro-benchmark, k=2 n=16): block solves of a fine sweep 3.2-3.8 ->
  1.1 ms, fine sweeps x1.11-1.24, level-1 sweeps x1.06-1.13 (~5% of the solve); results differ by ~1e-16 (relative).
  Using the inverse of the diagonal is allowed, but the LU solve also handles singular blocks (rank < block size),
  which an inverse doesn't: choose the method per block (inverse if `FullPivLU::rcond()` is large enough, LU
  otherwise), and check the iteration counts.

Other candidates:

- **Several aggregations at once (one P_F and one Galerkin product per level instead of one per pass)** — measured
  with existing options (2026-10-03, sequential): rejected as a default, no gain at k=0. The default aggregates can't
  be kept: the passes after the first aggregate on `A_T_Fc = Q_T^T A_T_F P`, where P_F depends on S (block Jacobi),
  so the intermediate Galerkin products are what the next pass sees, and they carry the anisotropy. Iterations (Cube-tet
  k=0 n=32 / Cube-cart-aniso100 n=32 / n=64 / Cube-tet-aniso20 n=16 / Cube-tet k=2 n=16), default 28 / 9 / 9 / 62 / 35:
  - aggregate with all neighbours at once, ignoring the strength (`-cs mn`): 35 / 39 / 49 / 60 / -;
  - passes on the cheap operators of Q_F, then one P_F (`-prolong 6 -coarsening-prolong 3`): 31 / 17 / 19 / 63 / 41.
    The averaged Q_F spreads the strong coupling of the removed faces over all the coarse faces: from the second pass
    on, the aggregates no longer follow the anisotropy (aniso100 n=32, level 1: 4109 aggregates and 23607 faces instead
    of 4096 and 23040). With `-face-prolong 2` (0 on the removed faces): 64 / 28 / - / 67 / -, but the same Q_F then
    also breaks the final reconstruction (`Theta()` no longer preserves the constants): to separate the two, the
    aggregation would need its own coupling operator (code);
  - passes driven by P_F, then one P_F (`-prolong 6`): 30 / 10 / 10 / 65 / 38: the one-shot P_F alone costs 1-3
    iterations;
  - `-cs dpa` (2 passes per level): 24 / 9 / 9 / 54 / -, but at k=0 7 levels and an operator complexity of 2.79 instead
    of 1.73 (solve 0.71 -> 0.96 s).

  Setup of a dedicated implementation of "passes on cheap operators, one P_F" (the setup of `-prolong 6
  -coarsening-prolong 3` minus the products of the passes that it computes uselessly, from the profile): k=0 ~3.6 s
  vs 3.55 s, no gain (the final P_F and Galerkin product on the fine level take 1.05 s: the one-shot P has as many
  non-zeros as the chained one, 2.10M vs 2.11M); k=2 ~6.3 s vs 8.6 s, but 41 iterations instead of 35 (setup + solve
  12.5 -> ~10.7 s). At most an opt-in for k>=1, if a cheap coupling operator that keeps the anisotropy is found.

- **Renumbering the fine unknowns for locality** (e.g. reverse Cuthill-McKee on the face graph) — measured on the
  kernels only (micro-benchmark, RCM): k=0 n=32, mean |i-j| 38000 -> 1900, fine SpMV 2.7 -> 2.0 ms, Gauss-Seidel
  sweeps -5 to -15%; k=2 n=16: ~-10%. Changes the order of the Gauss-Seidel sweeps, hence the iteration counts:
  opt-in only, and the gain is modest. Not worth it on its own.
- **Restriction without R** — measured (micro-benchmark, 2026-10-03): slower, saves memory only. `rc = P^T r` computed
  from the rows of P (row i adds `r_i P(i,:)` into rc, a scatter): R (= P^T, as many non-zeros as P) is then neither
  stored (25 MB at k=0 n=32, 115 MB at k=2 n=16) nor transposed in the setup (`R = P^T`: 0.24 s at k=2, see idea 5).
  But the scatter is slower than Eigen's `R * r`: k=0, 2.9 vs 2.0 ms sequential, 2.0 vs 0.9 ms on 16 threads
  (per-thread coarse vectors, then summed); k=2, 6.6 vs 6.3 ms and 4.9 vs 3.4 ms. Consecutive fine faces add into the
  same coarse entries (a dependency through memory), and the parallel version pays for its buffers and their sum.
  Keep R, unless memory runs short. This is not matrix-free: P is still stored (see the next item).
- **Matrix-free prolongation** (P not stored, applied from its factors) — estimated, medium effort, small gain.
  Structure of P_F (`-prolong 6`, `BuildProlongation`): P = E_K Q_F + E_R J Y, with Y = E_R Π Θ + E_K Q_F. E_K and
  E_R select the rows of the kept faces (41-43% of the fine faces; Q_F injects the coarse face: 1 non-zero per row)
  and of the removed faces (interior to an aggregate). Θ = -A_TcTc^-1 A_TcFc is the reconstruction in the coarse
  cells, and Π copies the constant mode of the coarse cell to the first coefficient of the removed face (only row 0
  of Θ is used). J = I - (2/3) D^-1 S is block Jacobi, used on the removed rows only. The product J Y fills in: a
  removed row of P has 11 non-zeros at k=0 and 62 at k=2, against 7 and 40 in the same row of S. Applying P e from
  the factors: Y (row 0 of Θ per coarse cell, injection on the kept faces), then `Y - (2/3) D^-1 (S Y)` on the
  removed rows. That reads ~25-35% fewer bytes than P (S's removed rows + the D^-1 blocks: 7 MB at k=2), but gathers
  in the fine vector Y, with the poor locality of the fine numbering, whereas `P * e` gathers in the coarse vector,
  which stays in cache. Estimate: prolongation -20-30% at k=2, perhaps more at k=0, where `P * e` (3.0 ms, 10 GB/s)
  is slowed by its many 1-non-zero rows. That is 2-4% of the solve. Applying P^T matrix-free is worse: a full product
  by S, plus scatters into the coarse vector (slower, see above). So restriction would keep `R * r`, and R is as
  large as P. Drawbacks: one implementation per prolongation option; the coarse mesh data (aggregates, Θ, D^-1) must
  survive the setup in compact arrays; the setup still builds P for the Galerkin product P^T S P (the peak memory
  doesn't drop). Rounding-level differences (the factors are applied in another order, and the coefficients <= 1e-15
  of P are no longer dropped): check the iteration counts. A cheaper first step, with the same gain on the 1-non-zero
  rows: store the kept faces of P as an index array (coarse face of each kept face: a copy, bit-identical) and only
  the removed rows in CSR.
- **Prolongation in place**: `x.noalias() += P*e` instead of `x += P*e` (a temporary): bit-identical, -5..10% of the
  prolongation (micro-benchmark), ~1% of the solve; needs a `Level::AddProlongation()` (U-AMG overrides `Prolong()`).

- `Level::ComputeGalerkinOperator()` computes `R * A * P` with Eigen's sequential products: used by the multigrids
  with the Galerkin operator that don't set their coarse operators themselves (e.g. the geometric multigrid with
  `-g 1`). `SparseMatrixOps::Multiply(SparseMatrixOps::Multiply(R, A), P)` gives the same matrix, bit for bit.
  Not measured yet: profile those setups first.
- C-AMG (`AggregAMG::GalerkinOperator`): sequential loop and triplets. Making it faster changes the U-AMG/C-AMG
  comparison of the paper (whose timings are tied to release 1.0 anyway).
- Memory: peak 1.33 GB for Cube-tet k=0 n=32 (whole run). The paper's n=64 extrapolates to 10-11 GB, too close to
  13 GB. Running the paper sizes here would first need a memory profile (probably the mesh and the assembly rather
  than the AMG: unverified).
