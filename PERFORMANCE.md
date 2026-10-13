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

AMD Ryzen 7 7735HS (8 cores, 16 threads), 13 GB of RAM, WSL2, a laptop: on battery, everything is ~35% slower
(2026-10-13); `cat /sys/class/power_supply/AC1/online` gives 1 when it is plugged in. No `perf` in WSL2; `valgrind` is installed (system,
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
  1. `ctest` (198 tests).
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
  the chunks) give results that depend on the number of threads, not on the run (since 2026-10-11: the values of the
  other chunks from the start of the sweep, see Done; before, a thread read x while another one wrote it). On 1 thread
  (and the scalar ones on a matrix of at most 20000 non-zeros), they are the sequential smoothers, bit for bit.

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

- **Library in Disk++** (2026-10-09/10, measurement only; Disk++ patch uncommitted, backup and scripts in
  `build/bench_lib/diskpp/`). Disk++'s `poisson_solver` and `poisson_tutorial` (u = exp(x) sin(pi x) sin(pi y)
  [sin(pi z)], K = Id; the tutorial: sin^3), tolerance 1e-10 (set by the clients), against its default direct solver
  (sequential MUMPS 5.8.2 + OpenBLAS; its time on 16 BLAS threads varies 2x from run to run: best of 16 threads / 1
  thread kept). Solver wall time (the client's: input conversion + setup + solve), fhhos4 on 8 threads, residual
  replacement every 20 iterations: 3D, tetrahedra 143k DoFs k=1 0.91 vs 1.73 s, 286k k=2 2.7 vs 6.3 s (60 its), 387k
  k=0 3.0 vs 5.0 s; 1.16M k=1 10.1 s, MUMPS out of memory (> 7.5 GB); Cartesian 286k k=1 1.37 vs 7.85 s, 774k k=0 6.2
  vs 20.5 s (1.8 vs 4.9 GB), 571k k=2 3.7 s, MUMPS out of memory. 2D: 0.8-1.2x on triangles (up to 1.57M DoFs),
  1.0-1.8x on quadrangles, 1.4-2.2x on hexagons (MUMPS is efficient in 2D). At k=0 the setup is 3-4x the solve (2D triangles 785k DoFs: setup 3.6 s, solve 1.0 s at 1e-8;
  Cartesian 3D 774k: 4.5 / 1.7 s). Iterations a bit above fhhos4's own meshes (Disk++ structured tetrahedra, 47k
  faces, 1e-8: 34/36/49 at k=0/1/2 vs 26/33/33 on the GMSH cube n=16), same time per iteration (k=1: 0.80 s / 36 its
  vs 0.83 s / 33 its in the program). With 1e-8: 3D 1.7-6.3x. The setup time printed at `Verbosity = 1` excludes the
  conversion of the inputs (column-major: ~0.07 s for 4.3M non-zeros).
  L2 errors vs MUMPS: <= 1e-6 relative everywhere (3D <= 3e-10), except 2D k=2 on fine meshes (L2 errors 3e-10 to
  3e-11): up to 1.2% (quadrangles 8 k=2), at any tolerance. Cause: Disk++'s condensed matrix is symmetric only up to
  rounding (|A - A^T|/|A| ~3e-16: local Schur complements) and fhhos4's FCG applies the lower triangle; MUMPS on the
  lower-symmetrized matrix gives fhhos4's L2 errors (quadrangles 7 k=2: 7e-8). Against it, 1e-8 leaves up to 2.7e-4,
  1e-10 1e-7..1e-6, 1e-12 no better (attainable true residual ~1e-11 in 2D, MUMPS 8e-12..3e-11). Logs `tol*/`.

- **Library: residual replacement in the FCG** (2026-10-10, the user's request). `FlexibleConjugateGradient::
  ResidualRecomputation` (0 in the program: its path unchanged; 20 in the library): the recursive residual is replaced
  by `b - Ax` every N iterations and when the stopping criteria are met (`IterativeSolver::StoppingCriteriaMet`, the
  criteria without the callback `OnIterationEnd`; `IterationResult::ReplaceResidualNorm`), and the loop condition
  evaluates them again on it: the iterations continue if they are no longer met. Near the attainable accuracy, the
  recurrence drifted from the true residual (Disk++ quadrangles 7 k=2, tolerance 1e-12: 8e-13 reported, 1.4e-11 true),
  so `Result::RelativeResidual` is now the true residual; a tolerance below the attainable accuracy runs
  `MaxIterations` (program's `-s fcglibuamg`, square cart n=64 k=2: 1e-13 met at 19 its, 1e-15 never, true residual
  ~2e-14). It did not change the L2 errors of the Disk++ runs (see above: the gap was the asymmetry). Cost: one
  matrix-vector product per 20 iterations and at convergence. `ctest` 132/132.
- **GCR and U-AMG for non-symmetric matrices** (2026-10-10, measurement only). GCR(10) (`-s gcr*`, library
  `Krylov = "gcr"`; its K-cycle iterates with GCR) on the symmetric matrices: 0-3 more iterations than FCG (Cube-tet
  n=8 k=1: 32 -> 35; square tri n=64 k=0: 20 -> 22). `SetupNonsymmetric()` (as AGMG: prolongations from the symmetric
  part, Galerkin `P^T A P`, LU on the coarsest level) on HArDCore3D's convection-diffusion (`hho-diff-advec-reac
  --solver_type fhhos4`, uncommitted in the fork), 1 thread, tolerance 1e-12, L2 errors equal to BiCGSTAB's
  everywhere (logs `HArDCore3D-release/build/runs_advec/advec_{16,32}.txt`; "solve" = HHO3D's solving time, with the
  block assembly and the setup):

  | Case (`-c`) | k | 16^3: BiCGSTAB / U-AMG its | 32^3: BiCGSTAB / U-AMG its | 32^3 solve: BiCGSTAB / U-AMG |
  |---|---|---|---|---|
  | `1 2 3 0` (anisotropic diffusion, advection (yz,xz,xy)) | 0 | 79 / 24 | 181 / 30 | 0.38 / 0.91 s |
  | | 1 | 137 / 19 | 259 / 23 | 4.19 / 2.93 s |
  | `1 1 5 0` (advection (1e5,0,0): dominated) | 0 | 32 / 26 | 64 / 27 | 0.16 / 0.73 s |
  | | 1 | 32 / 10 | 64 / 19 | 1.16 / 2.44 s |

  U-AMG: GCR(10) preconditioned by one K-cycle (GCR on the coarse levels), tolerance 1e-12 on the true residual.
  BiCGSTAB: Eigen's, diagonal preconditioner, default tolerance (machine epsilon): not at equal tolerance, which favors
  U-AMG. 32^3 k=1 `1 2 3 0`: U-AMG's 2.93-3.11 s = block assembly in `HHO3D::solve()` (sequential, ~0.55 s) + setup
  1.09 s + solve 1.47 s; HArDCore3D's assembly: 11.8 s (1 thread). Symmetric reference (HHO_Diffusion `-c 1 2`, same
  diffusion, 32^3 k=1, FCG): 13 iterations. Dominated case: h-robust at k=0 (26 -> 27), not at k=1 (10 -> 19).
  **Diagnosis of the dominated case** (2026-10-10; systems exported from HArDCore3D, scratch driver
  `build/bench_nonsym/driver.cpp` building U-AMG as the library does, logs `build/bench_nonsym/exp*.txt`): the
  ordering of the faces for the Gauss-Seidel smoother. In HArDCore3D's numbering (|strict upper|/|strict lower| of the
  condensed matrix 0.36-0.55), a sweep carries the transport only a few cells downstream; numbered along the flow
  (sorted by x: nearly triangular, 3e-4, or 3300 decreasing: the backward post-smoothing sweep is then the downwind
  one), U-AMG converges in 3-4 iterations at k=0 and k=1, 16^3 and 32^3 (32^3 k=1 solve 1.15 -> 0.23 s; setup 0.96 s).
  Two levels (k=1 + exact solve of the p-level k=0): 7 -> 10 iterations in HArDCore3D's order, 4 -> 4 along x: the
  k=1 modes only rely on the fine smoother, whereas at k=0 the first level is an aggregation along the flow, which
  compensates (26-27, flat). Restarts don't matter (GCR without restart: 19 -> 17); more smoothing halves the
  iterations without removing the growth (nu=2: 5 -> 9). Only the K-cycle converges in this case (V, W: no
  convergence, with GCR or BiCGSTAB); the symmetric part alone (FCG) is harder than the system (50 -> 79: strongly
  anisotropic, the upwind stabilization only on the faces normal to x). The diffusion-dominated case behaves like the
  symmetric one (2 levels: 16 -> 17; V-cycle converges).
- **Non-symmetric U-AMG: the matrix itself, averaged couplings, downwind order** (2026-10-11, the user's requests).
  No symmetric part is formed: the same algorithm runs on the matrix (`UncondensedLevel::Symmetric` false, `A_F_T`
  carried through the coarsening), without the optimizations of the symmetric case (lower triangular parts, Cholesky;
  first with global Galerkin products, see the next entry): reconstruction `-A_TT^{-1} A_TF` (inverse of the static
  condensation, LU), prolongation smoothing and Galerkin products with the matrix; couplings of the aggregation from the
  averaged traces (tr A_TF(T,F) + tr A_FT(F,T))/2 (equal to the symmetric ones for a symmetric matrix); Gauss-Seidel
  sweeps in the downwind order of each level (`Solver/Multigrid/DownwindOrder.h`; `SweepOrder` of `GaussSeidel`,
  `BlockGaussSeidel`), pre- and post-smoothing. The non-symmetric path on a symmetric matrix gives exactly the
  iterations of the symmetric one. HArDCore3D systems, scratch driver, 1 thread, iterations 16^3 -> 32^3
  (`build/bench_nonsym/exp11_nosym.txt`, `exp13_avg.txt`):

  | Case | first version (symmetric part) | row couplings (Notay), no order / + order | averaged couplings, no order / + order |
  |---|---|---|---|
  | dominated k=0 | 26 -> 27 | 38 -> 39 / 3 -> 4 | 19 -> 30 / 4 -> 4 |
  | dominated k=1 | 10 -> 19 | 14 -> 20 / 3 -> 4 | 12 -> 14 / 3 -> 4 |
  | diffusion-dominated k=0 | 24 -> 30 | 26 -> 36 / 26 -> 36 | 24 -> 29 / 24 -> 29 |
  | diffusion-dominated k=1 | 19 -> 23 | 18 -> 23 / 18 -> 23 | 18 -> 23 / 18 -> 23 |

  The row couplings (as Notay's strong couplings, from A_TF(T1,F) and A_FT(F,T2)) cost 20% at 32^3 in the
  diffusion-dominated case at k=0 (asymmetry |A - A^T|/|A| = 0.14%), already with 2 levels (37 vs 30): the averaged
  ones were kept (the user's choice). Setup cheaper than the first version (32^3 k=1: ~0.6 s instead of 1.05-1.45).
  Cost of the order: an ordered sweep is slower (irregular accesses, no residual deduced from the sweep: one more
  product per smoothing): ~0.09 s per iteration instead of ~0.06 (32^3 k=1); the topological sorts are sequential.
  Hybrid smoother (chunks of the order per thread): 4 iterations on 4 and 8 threads (17 without the order). Next: a
  cheaper topological sort; permuting the levels' matrices in the downwind order (locality, and the sweep kernels that
  deduce the residual).
- **Non-symmetric U-AMG setup: local Galerkin products and cheaper downwind order** (2026-10-11, the user's request:
  "the setup is already costly enough"). (1) `LocalMatrices` for a non-symmetric operator: `Decompose()` reads the full
  rows, `GalerkinProduct()` computes the full local matrices (no mirroring); exact (`tests/LocalMatricesTest.cpp`: equal
  to the global `P^T S P` to 1e-12 on an operator with 50% asymmetry). (2) `DownwindOrder()`: couplings in CSR arrays
  (computed once, no vector per row), upstream lists filled per thread by contiguous ranges of rows; same rule (fewest
  unprocessed upstream faces, then smallest number), hence the same order. A breadth-first order (FIFO, 2-7x faster) was
  tried and rejected: it puts all the faces without upstream ones first, and the hybrid smoother (chunks of the order
  per thread) then took 31-40 iterations instead of 4 on 8 threads. Fine level of the 32^3 systems, 1 / 8 threads:
  convection-dominated k=1 51 -> 41 ms / 42 -> 32 ms, k=0 41 -> 26 / 33 -> 24 ms, diffusion-dominated (no one-sided
  coupling) 32 -> 24 / 14 -> 12 ms. Setup before (global products, previous order) -> after, minimum of 3 alternated
  runs, library default smoothers (`build/bench_nonsym/bench_setup.sh`, `exp16_setup_before_after.txt`):

  | 32^3 | 1 thread | 8 threads |
  |---|---|---|
  | diffusion-dominated k=0 | 0.579 -> 0.520 s (-10%) | 0.373 -> 0.323 s (-13%) |
  | diffusion-dominated k=1 | 0.727 -> 0.620 s (-15%) | 0.457 -> 0.394 s (-14%) |
  | convection-dominated k=0 | 0.466 -> 0.535 s (+15%) | 0.348 -> 0.376 s (+8%) |
  | convection-dominated k=1 | 0.610 -> 0.658 s (+8%) | 0.443 -> 0.458 s (+3%) |

  Iterations unchanged in the diffusion-dominated case; in the convection-dominated one, the aggregation changed
  (without the order: 6 -> 5 levels at k=1, 30 -> 22 iterations at k=0), 3-5 iterations with the order either way.
  Not the rounding of the products, as first thought, but a bug of this change, fixed on 2026-10-12 (see "Bug of the
  non-symmetric path"): the coarse `A_F_T` was empty. Same binary, local / global products, time of the coarsening passes
  (prolongation + products) and setup, 1 thread: diffusion-dominated k=1 218 / 336 ms, 0.62 / 0.86 s;
  convection-dominated k=1 221 / 217 ms, 0.66 / 0.62 s; k=0 221 / 204 ms, 0.51 / 0.43 s: in the convection-dominated
  case the local products gain nothing in the passes (its coarse cells have slightly more faces: 4.43 vs 4.20 per cell
  on the first coarse level), and the local path also assembles the operator of each level from the local matrices.
- **Downwind order across the Péclet numbers** (2026-10-11, measurement only). The order is computed for every
  non-symmetric matrix (never with `SetupSymmetric()`), on each level but the coarsest, and kept where some couplings
  are one-sided: nothing to set by the user. Diffusion-dominated: it finds nothing (natural order), at the cost of the
  couplings' pass, 24-25 ms on 1 thread / 12-13 ms on 8 threads on the fine level (32^3, k=1). Sweep without a new test
  case: HArDCore3D's blocks `1 1 0 0` (diffusion Id) + gamma (`1 1 5 0` - `1 1 0 0`), i.e. the convection
  gamma (1e5, 0, 0), recondensed by the scratch driver (`base=`, `gamma=`; `build/bench_nonsym/sweep_peclet.sh`,
  `exp18_peclet.txt`), right-hand side A*1, 1 thread. 32^3, cell Péclet number gamma 3125; iterations and solve time,
  order off / automatic (ordered levels):

  | cell Péclet | k=0 | k=1 |
  |---|---|---|
  | 0.03-0.09 | 25-26 its / same (no level ordered) | 17 its / same (no level ordered) |
  | 0.31 | 26 its 0.20 s / 28 its 0.25 s (3/3) | 19 its / same (0/4) |
  | 0.94 | 27 its 0.20 s / 25 its 0.32 s (3/3) | 22 its 1.48 s / 22 its 1.38 s (2/4) |
  | 3.1 | 26 its 0.17 s / 23 its 0.24 s | 19 its 1.30 s / 17 its 2.14 s (4/4) |
  | 9.4 | 23 its 0.14 s / 17 its 0.19 s | 17 its 1.05 s / 13 its 1.39 s |
  | 31 | 24 its 0.14 s / 10 its 0.10 s | 18 its 1.09 s / 8 its 0.83 s |
  | 3125 | 24 its 0.19 s / 4 its 0.05 s | 18 its 1.10 s / 4 its 0.42 s |

  In iterations the automatic order is right everywhere (equal or fewer from a cell Péclet number ~1, at most 2 more
  around 0.3-0.6); in time it only pays from ~15-30, because an ordered sweep costs 1.7-1.9x a natural one (irregular
  accesses; on the scalar levels, no residual deduced from the sweep: one more product per smoothing). See idea 7.

- **Deterministic hybrid Gauss-Seidel** (2026-10-12, idea 6, the user's request). The hybrid smoothers (`hbgs`,
  `hrbgs`: scalar and block, natural order and `SweepOrder`; the library's default on several threads) are now
  Gauss-Seidel in the chunk of rows of each thread and a true Jacobi between the chunks: a row reads the current values
  of its chunk, and those of the other chunks from the start of the sweep (`src/Solver/FixedPoint/HybridChunks.h`).
  Before, a thread read the other chunks while they were written: the results varied from run to run, which cost much
  with GCR as the outer iteration and in the K-cycle: HArDCore3D's diffusion-dominated convection-diffusion
  (`-c 1 2 3 0`), 48^3 k=1, 8 threads, 191 iterations instead of 26 (sequential); 48^3 k=0, 105 instead of 46 (only
  GCR/GCR was that sensitive; on the symmetric part, FCG/FCG 20 -> 20). A prototype (the other chunks read in a copy
  of x, a test per non-zero) had given 24 iterations in every run at 32^3 k=1 on 4 threads (racy: 26-36).
  - Implementation: the chunks are contiguous ranges of positions of the sweep sequence (with `SweepOrder`, the chunk
    of a column through the order). Only the boundary block rows (a non-zero in another chunk) need old values: their
    non-zeros are described by ranges of positions in the row (external; internal before and after the diagonal
    block: in the natural order, a chunk is a range of columns, so at most two external ranges and one internal range
    on each side), found once per matrix, number of threads and order, at the setup (again if the number of threads
    changes); in the natural order, a block row is interior if the first and last columns of its rows are in its
    chunk. Before a sweep from x != 0, the threads copy the values of x of the exported block rows (those read by
    another chunk; split evenly) into a vector with x's indexing, then a barrier; the boundary rows multiply their
    external non-zeros by these old values, their internal ones by x. The interior rows keep the kernels of the
    sequential sweeps. From x = 0: no copy, no barrier, and the boundary rows skip the products by the values not
    computed yet, as the interior ones. On 1 thread, the sequential smoothers, bit for bit (the scalar one then also
    deduces the residual from the sweep: the "1 thread" item of idea 6).
  - Fraction of boundary rows / of external non-zeros, 4 chunks (8, 16): HArDCore3D's 32^3 fine level 13% / 3.2%
    (31% / 7.3%, 65% / 15.7%); U-AMG's levels of ~17000 faces 31-34% / 7.5-8.4%, of ~3500 faces 70-78% / 19-21%;
    fhhos4's Cube-tet n=32 (GMSH numbering) fine level 83% / 32% (92% / 43%), level 1 32% / 7.4%, level 2 63% / 15%.
  - Poorly separated chunks: on Cube-tet, the Jacobi between the chunks diverges (FCG + U-AMG, `-smoothers
    hbgs,hrbgs -threads 4`: k=0 n=32 and k=1 n=16 never converge; racy: 30-31 / 37-38 iterations; `bgs`: 28 / 33).
    A level whose external non-zeros exceed 20% of its non-zeros (`HybridChunks::MaxExternalFraction`) runs the
    sequential sweeps (the user's choice among three options): Cube-tet 28 / 33 iterations, as `bgs`. On 4 threads, it
    applies to the levels of a few thousand faces (19-21%: HArDCore3D's convection-dominated k=0, 23 -> 22
    iterations), never to the large levels of HArDCore3D's meshes. Rejected: the l1 hybrid Gauss-Seidel (Baker,
    Falgout, Kolev, Yang, SIAM J. Sci. Comput. 33 (2011): the diagonal of each boundary row plus the sum of |a_ij| over
    its external non-zeros; convergent for any SPD matrix and chunks): Cube-tet 35 / 72 iterations, HArDCore3D
    diffusion-dominated +1, convection-dominated 15 -> 19 (natural order), 5 -> 10 (downwind order).
  - Development (micro-benchmark of one sweep, HArDCore3D's 32^3 fine level, 4 threads, racy -> deterministic): a first
    version summed the external products of each chunk's boundary rows before the barrier, row by row: backward sweep
    k=1 6.6-7.6 -> 8.8 ms (sums 0.5-0.9 ms, gathering far values of x; the chunks between two others have twice the
    boundary rows of those at the ends: 0.45 ms of waiting at the barrier), forward sweep from 0, k=0, 0.19 -> 0.25 ms
    (the boundary rows multiplied the zeros not computed yet). The copy of the exported values instead: copy 0.005 /
    0.015 ms (k=0 / k=1), barrier 0.003-0.01 ms, the sweep itself as fast as the racy one (k=0 0.31-0.32 vs 0.32-0.35
    ms, k=1 6.7-6.8 vs 6.6-7.1 ms). Finding the boundary rows at the first sweep cost 1.3-2 ms (k=0) / 7-13 ms (k=1)
    per smoother and fine level: moved to the setup.

  Scratch driver (U-AMG built as the library does, GCR + K-cycle with GCR, HArDCore3D's 32^3 systems, smoothers
  `hbgs,hrbgs`; `build/bench_det/bench_driver.sh`, `driver_runs_final2.txt`; prototype: `driver_runs_final.txt`),
  iterations and solve time, racy / prototype / deterministic, 3 runs each (the prototype in a previous series; 3
  threads in the logs):

  | Case | 1 thread (`bgs,rbgs`) | 2 threads | 4 threads |
  |---|---|---|---|
  | diffusion-dominated k=0 | 29, 0.23 s | 31 / 31 / 31, 0.18-0.19 / 0.23-0.25 / 0.19-0.21 s | 30-33 / 31 / 31, 0.17-0.19 / 0.19-0.20 / 0.17-0.19 s |
  | diffusion-dominated k=1 | 23, 1.40 s | 23-24 / 23 / 23, 0.96-1.02 / 1.32-1.35 / 0.94-1.01 s | 30-31 / 24 / 24, 0.98-1.04 / 1.02-1.12 / 0.81-0.86 s |
  | convection-dominated k=0, natural order | 22, 0.16 s | 22 / 22 / 22, 0.13-0.14 / 0.15-0.16 / 0.14-0.15 s | 23 / 23 / 22, 0.11-0.13 / 0.14-0.15 / 0.11-0.12 s |
  | convection-dominated k=1, natural order | 15, 0.90 s | 15 / 15 / 15, 0.59-0.63 / 0.83-0.94 / 0.54-0.61 s | 15 / 15 / 15, 0.47-0.50 / 0.57-0.65 / 0.45-0.56 s |
  | convection-dominated k=0, downwind order | 4, 0.044 s | 5 / 5 / 5, 0.043-0.052 / 0.041-0.045 / 0.037-0.041 s | 5 / 5 / 5, 0.029-0.033 / 0.027-0.034 / 0.028-0.034 s |
  | convection-dominated k=1, downwind order | 4, 0.37 s | 4 / 4 / 5, 0.22-0.25 / 0.23-0.24 / 0.28-0.31 s | 4 / 4 / 5, 0.16-0.18 / 0.16-0.18 / 0.20-0.21 s |

  On 1 thread, the hybrid smoothers give exactly the results of `bgs,rbgs` (same iterations and residuals, 8 cases).
  The deterministic ones give the same residual in every run; the racy ones varied. With the downwind order, k=1, one
  more iteration on 2-4 threads (see idea 6). BiCGSTAB + V-cycle (diffusion-dominated, 4 threads): k=0 21 iterations
  (sequential: 21; racy: 21-25), k=1 28 (sequential 27; racy 26-31). HArDCore3D's case at 48^3 k=1 (`-c 1 2 3 0`,
  4 threads): 27 iterations, the same residual in 2 runs, solve 3.3-3.5 s (1 thread: 26 iterations, 5.8 s), setup
  1.8-2.1 s, L2 error unchanged, peak 4.79 GB. Symmetric matrices, 4 threads, FCG: HArDCore3D's `hho-diffusion -c 1 2`
  32^3, k=0 / 1 / 2: 30 / 19 / 26 iterations, as the racy version (3 runs each, the same residual in every run);
  fhhos4 (`-smoothers hbgs,hrbgs -threads 4`), Cube-cart k=1 n=32: 10 iterations, solve 0.25-0.29 s (racy: 10,
  0.25-0.26 s); Cube-tet k=0 n=32: 28, 0.58-0.74 s (racy 30-31, 0.57-0.58 s); Cube-tet k=1 n=16: 33, 0.62-0.69 s (racy
  38, 0.46-0.51 s: the fine level is sequential). `ctest` 168/168 (the program's default `bgs,rbgs` is unchanged),
  `quick_check.sh` OK.

  Against HArDCore3D's own solver (`hho-diff-advec-reac --solver_type bicgstab|fhhos4`: Eigen's BiCGSTAB, diagonal
  preconditioner, sequential, tolerance machine epsilon; fhhos4 1e-12), HArDCore3D's solving time (wall: for fhhos4,
  the blocks A_TT, A_TF, A_FT assembled from triplets in `HHO3D::solve()`, sequential, the conversion of these
  column-major inputs, the setup and the solve), one series (`build/bench_det/bench_speedup_det.sh`,
  `speedup_table_det.py`), L2 errors equal everywhere:

  | Case | BiCGSTAB | fhhos4 1 thread (setup + solve + other) | fhhos4 4 threads | speedup 1 / 4 threads |
  |---|---|---|---|---|
  | 32^3 k=1 diffusion-dominated | 211 its, 3.60 s | 23 its, 2.64 s (0.71 + 1.39 + 0.54) | 24 its, 1.90 s (0.53 + 0.82 + 0.55) | 1.36 / 1.89 |
  | 48^3 k=1 diffusion-dominated | 289 its, 17.5 s | 26 its, 10.9 s (3.09 + 5.71 + 2.10) | 27 its, 8.23 s (2.06 + 3.42 + 2.75) | 1.61 / 2.13 |
  | 32^3 k=1 convection-dominated | 61 its, 1.22 s | 4 its, 1.64 s (0.76 + 0.35 + 0.53) | 5 its, 1.35 s (0.56 + 0.21 + 0.58) | 0.74 / 0.90 |
  | 48^3 k=1 convection-dominated | 87 its, 5.92 s | 5 its, 7.37 s (2.78 + 1.72 + 2.87) | 6 its, 5.50 s (1.99 + 1.09 + 2.42) | 0.80 / 1.08 |
  | 32^3 k=0 diffusion / convection-dominated | 135 / 59 its, 0.29 / 0.17 s | 0.85 / 0.73 s | 0.74 / 0.48 s | 0.34 / 0.39, 0.23 / 0.35 |
  | 48^3 k=0 diffusion / convection-dominated | 189 / 84 its, 1.79 / 0.90 s | 3.72 / 2.60 s | 2.65 / 2.02 s | 0.48 / 0.68, 0.35 / 0.45 |

  From 1 to 4 threads, fhhos4's solve gets 1.6-1.7x faster at k=1, its setup 1.3-1.5x; the "other" part (the client's
  sequential triplet assembly and the conversion of the inputs: 2.1-2.9 s at 48^3 k=1) doesn't, and is a third to half
  of the time on 4 threads. Assembling the blocks in HArDCore3D's parallel assembly loop would give an estimated ~3x
  against BiCGSTAB at 48^3 k=1 diffusion-dominated, ~1.9x convection-dominated. At k=0, BiCGSTAB is cheaper than the
  setup alone.

- **Coarsening passes on flat arrays** (2026-10-12, idea 2 extended to the whole pass, the user's request: "Can we
  significantly improve setup time?"). `src/Solver/Multigrid/UncondensedAMG/CoarseningPass.h`, used with the default
  options (`UncondensedLevel::FlatPassesApply()`: `-face-prolong 1`, `-coarsening-prolong 6`, `-prolong 2`, `-cs
  mpa|dpa`, interface collapsing, local Galerkin products, no coarse `A_F_F` nor chained `Q_F`, LU diagonal blocks);
  the other options keep `HybridAlgebraicMesh`. Bit-identical.
  - A pass: the faces of the cells (rows of A_T_F), the cells of the faces (its transpose) and the neighbours sorted by
    coupling in CSR arrays (`IndexLists.h`), instead of vectors per element, face and aggregate, `std::map`s and
    mutexes; the same sequential pairwise aggregation on them; the coarse faces numbered in parallel by their pair
    (first aggregate, neighbour), the order of the interface collapsing. Then one parallel loop over the aggregates
    computes A_T_Tc, row 0 of Theta, the rows of P_F on the removed faces (the row of the operator merged from the
    local matrices of the two cells, its block Jacobi row, `J Y` in the order of the columns of J: idea 4, in the
    exact order of the global products) and A_T_Fc = Q_T^T A_T_F P: no Q_T, Q_F, Pi, triplets, `block()`,
    `IterationMatrix()` nor global products, and P is written by rows. The decomposition of the operator and the local
    Galerkin products are those of `LocalMatrices` (flat versions, which the mesh versions call). The 1x1 blocks (all
    the passes of `p_h`) reproduce the arithmetic of Eigen's 1x1 solves, which multiply by the inverse of the
    diagonal: LLT `-((b * (1/sqrt(a))) * (1/sqrt(a)))`, PartialPivLU `-(b * (1/a))`, block Jacobi
    `(s * (1/d)) * -omega`; the larger blocks call Eigen's (LLT, PartialPivLU, FullPivLU, reused per thread, and
    `BlockDiagonalSolver::SolveLU()`).
  - Checks: `uamg_harness compare` (operators and prolongations of all the levels, bit for bit) on Cube-tet k=0 n=32,
    k=2 n=16 with `-hp-cs h` (blocks 10/6) and `p_h`, 1 and 4 threads; HArDCore3D's 32^3 and 16^3 systems (scratch
    driver), symmetric (32^3, k=0, 1) and non-symmetric (diffusion- and convection-dominated, k=0, 1, with the downwind
    order): the same iterations and residual to 17 digits; `amg_papers.sh`: identical iteration tables (30
    runs); `ctest` 181/181, of which `tests/CoarseningPassTest.cpp` (13 cases: the program on 10 problems, sequential,
    the same iterations and L2 error bit for bit with `UncondensedLevel::UseFlatPasses = false`; the levels of a
    symmetric, rescaled and non-symmetric system).
  - Setup, before -> after, best of 2 (`build/bench_setup/`: `harness_flat` and `driver_flat`, the same binary on both
    paths, `FHHOS4_OLD_PASSES=1` for the old one; 4 threads: the library's `hbgs,hrbgs` in the driver):

    | Case | 1 thread | 4 threads |
    |---|---|---|
    | Cube-tet k=0 n=32 (318k faces) | 3.03 -> 1.29 s (2.4x) | 1.81 -> 0.61 s (3.0x) |
    | Cube-tet k=2 n=16, `p_h` (default) | 0.39 -> 0.23 s | 0.33 -> 0.12 s |
    | Cube-tet k=2 n=16, `-hp-cs h` | 3.99 -> 2.78 s | 3.09 -> 1.56 s |
    | HArDCore3D 32^3 k=0, symmetric | 0.51 -> 0.19 s (2.7x) | 0.37 -> 0.12 s (3.2x) |
    | HArDCore3D 32^3 k=1, symmetric | 0.62 -> 0.28 s | 0.45 -> 0.19 s |
    | HArDCore3D 32^3 k=0, non-symmetric (diffusion-dominated) | 0.54 -> 0.20 s | 0.39 -> 0.13 s |
    | HArDCore3D 32^3 k=1, non-symmetric (diffusion-dominated) | 0.65 -> 0.36 s | 0.47 -> 0.24 s |

    The program, its default threads (all logical cores, one per physical core in the solver), best of 2, before
    (`bench_setup/fhhos4_head`) -> after (`bench_setup/program_bench.sh`, logs in `program/`; the same iterations, the
    solve unchanged):

    | Case | Setup | Solve |
    |---|---|---|
    | Cube-tet k=0 n=32 (28 its) | 1.90 -> 0.52 s (3.7x) | 0.70 s |
    | Cube-cart-aniso100 n=64 (9 its) | 2.60 -> 0.53 s (4.9x) | 0.49 s |
    | Square-tri k=0 n=256 (24 its) | 1.33 -> 0.28 s (4.8x) | 0.37-0.44 s |
    | Cube-tet k=1 n=16, `p_h` (33 its) | 0.24 -> 0.09 s (2.7x) | 0.69 s |
    | Cube-tet k=2 n=16, `-hp-cs h` (35 its) | 2.09 -> 1.35 s (1.5x) | 3.5 s |

    At k=0 the setup is now shorter than the solve (HArDCore3D 32^3: 0.19 vs 0.24 s on 1 thread). At k=1, the
    extraction of the p-level and the smoothers of the k=1 level (~0.1 s on 1 thread) are now a third of the setup.
    The papers' runs (`amg_papers.sh`, default threads, one run each): Fig. 4.4 Cube-tet n=32 1.94 -> 0.56 s,
    Cube-cart n=32 0.31 -> 0.08 s, Complex-tet n=4 0.56 -> 0.14 s, Heterog1e8 n=256 0.37 -> 0.07 s.
  - What remains (Cube-tet k=0 n=32, 1 / 4 threads, timers of the prototype): the loop over the aggregates 0.39 /
    0.11 s (latency of the accesses to the data of the fine level, in the numbering of the mesh: gathering the data of
    an aggregate 0.13 s, merging the rows of the local matrices 0.13 s on 1 thread), chaining of P
    (`Multiply(P, P_pass)`) 0.22 / 0.12 s, local Galerkin products 0.19 / 0.06 s, neighbours and couplings 0.11 /
    0.05 s, interface collapsing 0.10 / 0.04 s, decomposition 0.08 / 0.03 s, pairwise aggregation (sequential) 0.05 s,
    cells of the faces (sequential transposition) 0.03 s.
  - How it was developed: a prototype of the whole pass, bit-identical at the first try (1.9x on 1 thread), then
    without allocations in the loop over the aggregates (Eigen objects and buffers reused per thread, the rows merged
    instead of sorted: 2.1x), then the 1x1 blocks without Eigen's dynamic-size solves (each a call to its blocked
    triangular solve, with its blocking space: 2.4x).
- **Bug of the non-symmetric path, fixed** (2026-10-12, found by the flat passes). With the local Galerkin products
  (2026-10-11, `bbb8e073`), `UncondensedLevel::Coarsen()` only transposed P for the coarse `A_F_F` or the global
  products, so the coarse `A_F_T = P^T A_F_T Q_T` was computed with an empty `Pt` (0 x 0): an empty matrix in Release
  (an assertion in Debug). The passes after the first and the coarser levels then computed the couplings with the
  traces of A_F_T as 0, i.e. from A_T_F alone, halved: the same aggregates for a symmetric matrix (`-uamg-nonsym`),
  hence unnoticed by the tests. Fixed in both paths (`Pt` computed for a non-symmetric matrix). HArDCore3D's 32^3
  systems, 1 thread: the same iterations with the downwind order (diffusion-dominated k=0 29, k=1 23;
  convection-dominated k=0 4, k=1 4, with 6 levels instead of 5); without the order, convection-dominated k=0 22 -> 30
  iterations, k=1 15 -> 15 (5 -> 6 levels): the "6 -> 5 levels, 30 -> 22 iterations" of the entry "Non-symmetric
  U-AMG setup" were this bug. The averaged couplings, chosen by the user on 2026-10-11, had been measured before it.

- **Flat passes for the other options** (2026-10-13, the user's request: "start moving the other options"). The flat
  passes (`UncondensedLevel::FlatPassesApply()`) now run with the face prolongations 1, 2 (`-face-prolong`), the
  coarsening prolongations 3 to 6 (`-coarsening-prolong`: all those of the local Galerkin products), the aggregations
  `-cs mpa|dpa|n|mn` (with all the neighbours: `CoarseningPass::Aggregate(true)`), the multigrid prolongations other
  than 7 (`-prolong 1|3|4|5|6|8`, besides the default chained ones: the aggregation composed over the passes, then
  `BuildProlongation()` on the fine mesh of the level, as before), the coarse `A_F_F` (p-levels after h-levels), the
  chained `Q_F` (coarse operators `Q_F^T A Q_F`, not Galerkin) and `-diag-block-solve inv`. `BuildProlongation()` and
  its traces now take a `FaceAggregation` (the removed faces, their aggregate, the cells of the faces), built from
  `HybridAlgebraicMesh` or from the flat passes: one implementation for both. Bit-identical (`uamg_harness compare`,
  `build/bench_setup/variant.sh`: every variant above on Cube-tet k=0 n=32 and k=2 n=16 `-hp-cs h`, and 24 small
  dumps of `build/local-galerkin/dump/check/`: polygonal, biharmonic, k=1-3 meshes...). Setup, old -> flat, 4 threads,
  one run each (indicative):

  | Cube-tet k=0 n=32 | old -> flat | Cube-tet k=2 n=16, `-hp-cs h` | old -> flat |
  |---|---|---|---|
  | `-coarsening-prolong 3 / 4 / 5` | 1.86 / 1.81 / 1.87 -> 0.53 / 0.63 / 0.58 s | `-coarsening-prolong 3 / 4 / 5` | 1.53 / 2.87 / 1.97 -> 0.64 / 1.89 / 0.64 s |
  | `-face-prolong 2` with 3 / 4 / 5 / 6 | 1.57 / 1.87 / 1.94 / 2.00 -> 0.41 / 0.59 / 0.51 / 0.61 s | `-diag-block-solve inv` with 6 / 4 | 3.08 / 2.77 -> 1.90 / 1.64 s |
  | `-cs dpa / n / mn` | 2.04 / 1.65 / 1.56 -> 0.56 / 0.54 / 0.62 s | `-cs mn` | 2.93 -> 1.94 s |
  | `-prolong 4 / 5 / 6 / 8` | 3.12 / 3.00 / 3.03 / 3.42 -> 1.37 / 1.30 / 1.25 / 2.08 s | `-prolong 6` | 5.84 -> 4.36 s |
  | coarse operators `Q_F^T A Q_F` | 2.61 -> 1.05 s | `-hp-cs h_p -num-meshes 2` (coarse `A_F_F`) | 7.38 -> 6.12 s |

  `amg_papers.sh`: identical iteration tables (30 runs); its runs of the variants of Fig. 4.5 (`-coarsening-prolong
  3|4|5`, Cube-tet n=16 and Cube-cart-aniso100 n=32, default threads): setup 0.20-0.33 -> 0.05-0.07 s. `ctest` 196/196
  (`CoarseningPassTest.cpp`: 28 cases).
  The non-chained prolongations and the coarse `A_F_F` keep their global products (`P^T A P` of the level, `P^T A_F_F
  P`): most of what remains of their setup. Still on `HybridAlgebraicMesh` only: the coarsening prolongations 1, 7, 8
  (global Galerkin products in the passes), `-prolong 7`, the agglomeration of the removed faces (`-fcs i`) and
  `-face-prolong 3`.
  Found in the old path while comparing (both paths now behave the same): `-prolong 3` crashed (`free(): invalid
  pointer`: `BuildProlongation()` returns the level's own `Q_F`, which was then moved into P and deleted; now copied,
  `SetProlongation()`); `-prolong 9` returned an uninitialized pointer (its body is commented out; now an error).
  `-prolong 1` makes the coarse blocks dense: a coarse face then belongs to more than two coarse cells, which the
  interface collapsing doesn't handle (the old path crashed on the second level: faces left without coarse face, and a
  data race on `CoarseElements` in its loop over the removed faces); the flat passes now stop with an error there. It
  was already impractical: 431 s of setup on Cube-tet k=0 n=32 (2 levels). `-face-prolong 3` segfaults in the old path
  (Cube-tet n=8), unchanged.

- **Default options: the rest of the setup** (2026-10-13, the user's request: "improve the performance of the default
  options"). Profile first (instrumented copies in `build/bench_setup2|3/`): after the flat passes, the setup of
  Cube-tet k=0 n=32 on 8 threads was spread over ~15 steps of 0.01-0.09 s, of which 0.21 s sequential (Eigen's
  transpose `R = P^T`, the pairwise aggregation, the cells of the faces); at k=1 (`p_h`, HArDCore3D 32^3, 8 threads,
  0.17 s), the level k=1 took half of it: the setups of its two smoothers (0.063 s: each factorized the same diagonal
  blocks and found the same hybrid chunks) and the extraction of the p-level (0.035 s, triplets). All bit-identical
  (levels compared with `uamg_harness compare`, solves with the scratch driver's residuals to 17 digits, 1 and 8
  threads, hybrid smoothers included):
  - `SparseMatrixOps::Transpose()` in parallel (from 20000 non-zeros: each thread counts, then copies, a contiguous range
    of rows: Eigen's result), used for `R = P^T` (`UncondensedLevel::SetupRestriction()`) and all the transposes of the
    setup. `R = P^T` on 8 threads: 0.083 -> 0.014 s.
  - The post-smoother reuses the setup of the pre-smoother when they are of the same kind on the same matrix
    (`Smoother::Setup(A, preSmoother)`, `IterativeSolver::SetupFrom()`: `BlockGaussSeidel` copies the factorized
    diagonal blocks, their positions and, with the same order, the hybrid chunks; `GaussSeidel` its diagonal and
    chunks). k=1 level, 8 threads: post-smoother 0.03 -> 0.003 s.
  - The p-level blocks extracted by rows (`ExtractCoarseMatrix()`: `BuildByRows`, no triplets): 0.035 -> 0.008 s.
  - The cells of the faces in parallel (atomic counts, then the cells of each face sorted: 2 at most): 0.025 -> 0.011 s
    on 8 threads (+0.005 s on 1 thread).
  - The pairwise aggregation (sequential: its order defines the aggregates) with the data of each cell packed together
    and prefetched a few cells ahead along the aggregation order: 0.053 -> 0.042 s (1 thread), 0.071 -> 0.053 s (8).
  - `Multiply()`: a row of A made of a single 1 (the injections of the prolongations) copies the row of B. Chaining of
    P on 8 threads: 0.092 -> 0.072 s.

  Setup, before (`0831fb25`) -> after, plugged in, medians of 5 interleaved runs (`harness_flat`, `driver_flat`; the
  step timings above were measured on battery, ~35% slower):

  | Case | 1 thread | 8 threads |
  |---|---|---|
  | Cube-tet k=0 n=32 | 1.17 -> 1.15 s (-2%) | 0.47 -> 0.38 s (-20%) |
  | HArDCore3D 32^3 k=0 | 0.173 -> 0.170 s (-2%) | 0.094 -> 0.081 s (-14%) |
  | HArDCore3D 32^3 k=1 | 0.27 -> 0.23 s (-14%) | 0.17 -> 0.11 s (-36%) |

  The program at its default threads, medians of 3, the same iterations (since `27ea0d44`, before the flat passes, in
  parentheses): Cube-tet k=0 n=32 0.47 -> 0.39 s (1.90 s), Cube-tet k=1 n=16 0.089 -> 0.061 s (0.24 s), Square-tri k=0
  n=256 0.24 -> 0.19 s (1.33 s), Cube-cart-aniso100 n=64 0.51 -> 0.46 s (2.60 s). `ctest` 198/198, `amg_papers.sh`:
  identical iteration tables.

  Sequential steps run slower when 8 threads are in use (Eigen's transpose alone: 0.046 s on 1 thread, 0.092 s on 8):
  after a parallel region the other threads spin a while (libgomp), and, on battery, the clock of the cores drops
  with their number. What remains (Cube-tet k=0 n=32, 8 threads): the loop over the aggregates 0.08 s, chaining of P
  0.07 s, interface collapsing, decomposition, local Galerkin products, neighbours and couplings 0.05-0.06 s each,
  the aggregation 0.05 s; on HArDCore3D's hexahedra, no step above 15% of the setup.

## Current profile (U-AMG; parallel: 16 threads, the default before 2026-10-08)

| Case | Setup, sequential / parallel | Solve, sequential / parallel |
|---|---|---|
| Cube-tet k=0 n=32 | 3.2 / 1.7 s | 0.7-0.9 / 0.7-0.9 s |
| Cube-cart-aniso100 n=64 | 4.2 / 2.5 s | 0.5 / 0.5 s |
| Cube-tet k=2 n=16 | 3.7 / 1.8 s | 3.7 / 3.5 s |

**Setup breakdown before the flat passes** (2026-10-12, see "Coarsening passes on flat arrays" in Done for what
remains; timers of `Prof.h` in a copy of `src/`, `build/bench_setup/`:
`uamg_harness` on the dump of Cube-tet k=0 n=32, the scratch driver on HArDCore3D's 32^3 systems; logs in `logs/`).
The setup is the h-coarsening passes of the first h-coarsened level, at any k: at k >= 1 (`p_h`), the p-levels only
cost their extraction and smoothers, and the passes run on the k=0 level.
- Cube-tet k=0 n=32 (318k faces), 1 / 4 threads: setup 3.57 / 1.94 s, of which level 0 3.30 / 1.77 s: 4 passes on
  nearly the fine size (faces 1.45x, 2.22x, 3.45x, 5.57x fewer; tetrahedra share one face per pair). By step, all
  levels: hybrid meshes 1.47 / 0.80 s (`Build()` 0.55 / 0.23, of which neighbours and couplings 0.36 / 0.12;
  `Coarsen()` 0.85 / 0.52, of which the aggregates' neighbours (`std::map`) 0.28 / 0.11, removed faces 0.21 / 0.08,
  pairwise aggregation 0.18 / 0.16 and interface collapsing 0.17 / 0.17, both sequential; auxiliary coarse meshes
  0.08 / 0.05); P_F 0.73 / 0.32 s (`IterationMatrix()` 0.18 / 0.05, rows of the removed faces and their block Jacobi
  0.17 / 0.07, `Theta()` 0.16 / 0.08, `J * Y` 0.09 / 0.04, `Pi`, `Pi * Theta`, row selections 0.12 / 0.06); local
  Galerkin products 0.51 / 0.17 s (products 0.32, decomposition 0.15, assembly 0.04 on 1 thread); `BuildQ_F`,
  `BuildQ_T` 0.24 / 0.15 s; the other products (`Q_T^T A_T_F`, `A_T_Tc`, `A_T_Fc`) 0.21 / 0.09 s; chaining of P
  0.14 / 0.08 s; freeing the meshes and matrices of the passes 0.10 / 0.15 s; `R = P^T`, smoothers, coarse solver
  0.04 / 0.04 s. Peak 0.6 GB.
- HArDCore3D 32^3, symmetric path (driver, `sym=1`, FCG), 1 / 4 threads: k=0 (95k faces) setup 0.53 / 0.36 s for a
  solve of 0.22 / 0.16 s, of which the 3 passes of level 0 0.45 s (1 thread); k=1 setup 0.71 / 0.45 s (solve 1.02 /
  0.67 s): extraction of the p-level 0.05 / 0.03 s, smoothers 0.06 / 0.04 s, passes of the k=0 level 0.51 / 0.29 s.
  Non-symmetric path (`diff_n32_k1`, 1 thread): 0.73 s, of which extraction 0.09 s, downwind orders 0.06 s, passes
  0.44 s. fhhos4 Cube-tet k=2 n=16 (`p_h`, 1 thread): 0.46 s, of which smoothers and extraction 0.14 s, passes 0.30 s.
- No step dominates (about 15 steps of 0.1-0.4 s on 1 thread), and none is arithmetic: a pass costs ~10 us per
  element at level 0 (2.9 s for the 305k elements of the 4 passes), thousands of cycles for 4-10 faces, i.e.
  vectors per element, face and aggregate, `std::map`s, mutexes, triplets sorted by `setFromTriplets`, `block()`
  scans and the freeing of all of it. On 4 threads, the sequential parts (pairwise aggregation, interface collapsing,
  freeing) take 0.48 of the 1.94 s. Hence idea 2, extended to the whole pass (done).

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
2. **Flat data structures for the hybrid algebraic meshes** (`HybridAlgebraicMesh`) — done for the default options
   (2026-10-12, "Coarsening passes on flat arrays" in Done: setup 2.4x faster on 1 thread, 3x on 4, at k=0). Next
   steps: the chaining of P and the locality of the loop over the aggregates (see what remains there). Most other
   options followed on 2026-10-13 ("Flat passes for the other options"); still on `HybridAlgebraicMesh`:
   `-coarsening-prolong 1|7|8`, `-prolong 7`, `-fcs i`, `-face-prolong 3`. The notes below are those of the estimate. Main lever at k=0 (the paper's case). Vectors of vectors, a mutex per element, a `std::map` of neighbours
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
   Re-measured 2026-10-12 (see "Setup breakdown" in Current profile), and extended to the whole pass: at k=0, the
   adjacency of the cells and faces is the pattern of A_T_F (CSR) and its transpose; the coarse faces numbered by a
   parallel sort of the pairs of aggregates (the order of the interface collapsing: by the first aggregate, then the
   neighbour); Q_T, Q_F, Pi and P_F (idea 4) written directly in CSR, without triplets; `Theta()` per aggregate (idea
   3); buffers reused from one pass to the next instead of freed. That covers ~75% of the setup (meshes 1.47, P_F
   0.73, Q_T/Q_F 0.24, products 0.21, freeing 0.10 of 3.57 s on 1 thread). Estimate, unverified: setup 2-3x faster on
   1 thread, more on 4 threads (the sequential pairwise aggregation over flat arrays: a few tens of ms). Large effort,
   in steps that can each be bit-identical (same aggregates, numbering and order of the sums).
3. **`Theta()` without `SparseMatrix::block()`** — done in the flat passes (default options); the other options still
   use `HybridAlgebraicMesh::Theta()`. Measured, small effort, bit-identical if the blocks and the solves are the same. `HybridAlgebraicMesh::Theta()` (the reconstruction of P_F, `-A_TcTc^-1 A_TcFc` per coarse cell)
   extracts its dense blocks with `block()`, which scans the row from its start for each face: 0.53 / 0.29 s of the
   3.8 / 1.8 s of setup at k=2 (sequential / 16 threads: it scales poorly), 0.15 / 0.07 s at k=0 (`uamg_harness`,
   2026-10-04). One pass over the rows of the cell, as in `BuildElementFaces()`, gives all the blocks; then a
   parallel fill without triplets. `CouplingValue()` (`-fcs i`) and `BuildHighOrderTraceOnRemovedFaces()` do the same.
4. **Local P_F** — done in the flat passes, bit-identical (the sums in the order of the global products). Estimated,
   medium effort, rounding-level differences. The rows of P_F on the removed faces only
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
   - Determinism, and the sequential sweeps on 1 thread: done (see Done, "Deterministic hybrid Gauss-Seidel").
   - In the downwind order, the transport crosses one chunk boundary per sweep (Jacobi between the chunks): 5
     iterations instead of 4 at k=1 on 2-4 threads (HArDCore3D's convection-dominated case, 32^3; the racy version
     sometimes read updated upstream values). Recovering it would need the threads to follow the flow from one chunk
     to the next (a wavefront: less parallelism); not tried.
7. **Renumbering the levels of a non-symmetric matrix in their downwind order** — measured in part, medium effort.
   Instead of sweeping in the downwind order (`SweepOrder`), permute each level that gets an order once at the setup
   (its operator, the rows of its prolongation and the columns of the finer one, the restrictions; the fine level's
   vectors at the entry and exit of the solve): the smoothers then sweep in their natural order, with the locality and
   the residual deduced from the sweep of the natural kernels, and `SweepOrder` goes away. Gain: an ordered iteration
   costs 0.09-0.1 s at 32^3 k=1 (convection-dominated), a natural one ~0.06 s (the fine faces sorted along x before
   the setup, scratch driver: 4 iterations at 0.056 s): ~0.04 s per iteration, and the order would pay in time from a
   cell Péclet number ~1 instead of ~15-30 (see Done, "Downwind order across the Péclet numbers"). Cost, measured on
   the 32^3 fine level, P A P^T: Eigen's generic permutation 470-500 ms at k=1 (18-32 ms at k=0): not usable; a
   hand-written one (rows copied to their new position, columns relabeled and sorted per row, parallel; checked:
   B x_new = (A x)_new to 1e-16) 135 ms on 1 thread / 28 ms on 8 threads at k=1, 9 / 3 ms at k=0. Sorting the columns
   once per block row (same pattern for the rows of a face) should bring k=1 to ~40-50 ms on 1 thread (copying the
   9.4M non-zeros is the floor; estimate), plus the coarse levels and the prolongations (~+20-30%). Only on the levels
   that get an order (never symmetric or diffusion-dominated). Not done (the user's decision, 2026-10-11: keep it as a
   possibility).

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

- **Coarsening prolongation without its block-Jacobi smoothing** (`-coarsening-prolong 5`, P_F^(0), instead of the
  default 6, P_F) — measured (2026-10-12), rejected (the user's decision: the anisotropic case worsens too much). It
  skips `BlockJacobi::IterationMatrix()` and the product `J * Y` in each pass, and P gets sparser; the operator
  complexity doesn't change (the pattern of the coarse operators comes from the aggregates). FCG + U-AMG, the program's
  defaults otherwise, best of 2 (`build/bench_prolong/`, `summary.txt`), iterations and setup / solve / total ratios
  P_F^(0) / P_F, 1 thread (4 threads):
  - Cube-cart k=0 n=32: 15 -> 15 its, setup 0.91 (0.89), total 0.92 (0.91); Cube-tet k=0 n=32: 28 -> 29, setup 0.89
    (0.90), solve 1.14, total 0.94 (0.96); Complex-tet n=4: 25 -> 25, total 0.88 (0.92); Cube-tet-aniso20 n=16: 62 ->
    63, total 0.93 (0.96); Square-tri k=0 n=512: 26 -> 23, total 0.91 (0.95); Square-tri k=1 n=256: 17 -> 17, total
    0.95 (0.94); Cube-tet k=1 n=16: 33 -> 33, total 1.03 (0.99); Cube-tet k=2 n=16 (4 threads): 33 -> 33, total 0.97;
    Cube-cart k=1 n=16: 9 -> 9, 0.08-0.09 s either way;
  - Cube-cart-aniso100 n=32: 9 -> 16 its, solve 1.98 (1.88), total 0.98 (1.05); Heterog1e8 n=256: 22 -> 27 its,
    solve 1.21, total 1.00 (1.12).
  Setup -5 to -15% at k=0 (at k >= 1, the passes run on the k=0 p-level only: -3 to -9%), total -4 to -12% on the
  isotropic homogeneous cases, but the iterations nearly double on the strong anisotropy and grow by 23% on the
  heterogeneous case, as in the paper's comparison of the variants.

- **MUMPS on a larger coarsest level** — measured (2026-10-12), rejected for now: no gain on HArDCore3D's 32^3
  systems. MUMPS 5.8.2 (sequential library, automatic ordering, OpenMP + OpenBLAS threads) as the coarse solver of the
  scratch driver (`MumpsSolver` in `build/bench_nonsym/driver.cpp`, `driver_mumps`, linked to the conda env `diskpp`'s
  MUMPS; analysis + factorization at the setup, a solve per visit; `SYM=1` on the symmetric path), against Eigen's
  `SparseLU` / `SimplicialLDLT`; same iterations and residuals to 3-4 digits. The coarsening divides the size by
  ~4-5 per level, so the coarsest-size limit (`MatrixMaxSizeForCoarsestLevel`, 1000) of 10000 / 30000 / 100000 removes
  one / two levels / all the multigrid at k=0 (at k=1, the coarsest is then the p-level k=0). Setup + solve, best of 2,
  4 threads (`build/bench_mumps/bench_coarse_mumps.sh`, `runs.txt`, `summarize_mumps.py`), relative to today's
  (Eigen, <= 1000 rows):

  | System | MUMPS, today's levels | MUMPS ~3200-3800 rows | MUMPS ~16600-18300 rows | MUMPS ~95000 rows | Eigen ~3200-3800 rows |
  |---|---|---|---|---|---|
  | diffusion-dominated k=0 (0.50 s) | 1.05 | 1.10 | 1.50 | 2.29 (direct solve) | 1.45 |
  | diffusion-dominated k=1 (1.30 s) | 0.99 | 1.03 | 1.19 | 1.67 (17 its instead of 24) | 1.17 |
  | convection-dominated k=0 (0.44 s) | 0.95 | 0.97 | 1.43 | 2.57 (direct solve) | 1.70 |
  | convection-dominated k=1 (0.80 s) | 0.89 | 0.94 | 1.16 | 2.05 | 1.15 |
  | symmetric part k=0, FCG (0.50 s) | 1.01 | 1.10 | 1.36 | 1.90 (direct solve) | 1.45 |
  | symmetric part k=1, FCG (1.18 s) | 0.94 | 0.91 | 1.08 | 1.55 | 1.10 |

  1 thread: the same picture (MUMPS ~3300 rows 0.96-1.13, ~17000 rows 1.15-1.41, direct 1.52-2.75). MUMPS factorizes
  10-35x faster than Eigen (~3500 rows: 10-23 vs 40-275 ms; ~17000 rows: 0.15-0.24 s vs 2.6-8 s in an earlier sweep)
  and solves 2-3x faster, but the levels a larger coarsest level removes are cheap (setup 0.01-0.04 s), whereas it
  costs its factorization and a solve at each visit: the K-cycle visits the coarsest level 30-180 times per solve in the
  diffusion-dominated cases, 0.8-1.2 ms per MUMPS solve at ~3300 rows, 3-5 ms at ~17000 (2-3x the time of solves back to
  back: the factors are evicted between visits). MUMPS on the whole fine system at k=0: its factorization alone (0.9-1.4
  s) exceeds U-AMG's setup + solve. It would matter if the coarsening stopped early on a large coarsest level, where
  Eigen takes seconds; with today's levels, MUMPS is within the noise of Eigen (0.89-1.05).

  In the program (a temporary `ENABLE_MUMPS` option, removed after the measurement: coarse solver `mumpsch`, LDL^T
  without pivoting; MUMPS 5.7.3 from conda-forge's `mumps-seq`, without OpenMP: its threads are OpenBLAS's,
  `OPENBLAS_NUM_THREADS` set to the number of threads; the input copied to MUMPS's coordinate format, 1-based, one
  triangle when symmetric: it takes no CSR), FCG +
  U-AMG with the default smoothers (`bgs,rbgs`), setup + solve relative to today's (`ch`, <= 1000 rows), one run
  (2D, 4 threads: best of 3) (`build/bench_mumps/bench_program_mumps.sh`, `program/`, `parse_program_mumps.py`):

  | Problem (today's setup + solve) | threads | MUMPS, today's levels | MUMPS, 1-2 levels less | MUMPS, 2-3 levels less |
  |---|---|---|---|---|
  | Cube-tet k=0 n=32 (4.71 / 3.10 s) | 1 / 4 | 1.08 / 1.07 | 9998 rows: 1.03-1.07 / 1.05-1.08 | 57111 rows: 1.53 / 1.84 |
  | Cube-tet k=1 n=16 (1.11 / 0.96 s) | 1 / 4 | 1.00 / 1.01 | 7666 rows: 1.05-1.08 / 1.14-1.25 | 43656 rows: 1.17 / 1.43 |
  | Square-tri k=0 n=512 (15.0 / 9.47 s) | 1 / 4 | 1.09 / 0.97 | 7169 rows: 1.02 / 0.93 | 28207 rows: 0.98 / 0.94 |
  | Square-tri k=1 n=256 (4.38 / 3.13 s) | 1 / 4 | 0.93 / 0.96 | 7045 rows: 0.99 / 0.93 | 51859 rows: 1.07 / 1.05 |

  The same iterations (+-1). In 3D, a larger coarsest level is slower; in 2D (where the fill-in of a sparse direct
  solver is low), 0.93-0.94 on 4 threads at 7000-28000 rows, of which 0.96-0.97 come from MUMPS at today's size: at
  the edge of the noise (runs within 3-10% of each other). Dropped (the user's decision).

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
