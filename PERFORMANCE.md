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

## Current profile (U-AMG, default OpenMP settings)

| Case | Setup, sequential / parallel | Solve, sequential / parallel |
|---|---|---|
| Cube-tet k=0 n=32 | 3.7 / 2.0 s | 0.8 / 0.8 s |
| Cube-cart-aniso100 n=64 | 4.7 / 2.8 s | 0.5 / 0.5 s |
| Cube-tet k=2 n=16 | 8.5 / 3.7 s | 3.7 / 3.6 s |

- **The solve barely benefits from the threads.** The sweeps of the default smoother (lexicographic Gauss-Seidel, as
  in the paper) are sequential; only the residual/`Ax` products and the intergrid transfers are parallel. Those are
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
   the smoothers, after which 16 threads are no longer slower than 1 at k=2 (3.6 vs 3.7 s): re-measure.
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
   Gauss-Seidel (e.g. multicolour) would be needed; its gain is unknown. Measured before the optimization of the
   smoothers, which made the sequential sweeps ~2x faster at k=2 (and the rows of the hybrid one too): re-measure.
   The hybrid smoother also has a data race (see the exactness facts): a real Jacobi between the chunks (or a
   colouring of the chunks' boundaries) would make it deterministic.

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

- **Renumbering the fine unknowns for locality** (e.g. reverse Cuthill-McKee on the face graph) — measured on the
  kernels only (micro-benchmark, RCM): k=0 n=32, mean |i-j| 38000 -> 1900, fine SpMV 2.7 -> 2.0 ms, Gauss-Seidel
  sweeps -5 to -15%; k=2 n=16: ~-10%. Changes the order of the Gauss-Seidel sweeps, hence the iteration counts:
  opt-in only, and the gain is modest. Not worth it on its own.
- **Restriction without R** — measured (micro-benchmark, 2026-10-03): slower, saves memory only. `rc = P^T r` computed
  from the rows of P (row i adds `r_i P(i,:)` into rc, a scatter): R (= P^T, as many non-zeros as P) is then neither
  stored (25 MB at k=0 n=32, 115 MB at k=2 n=16) nor transposed in the setup (`R = P^T`: 0.24 s at k=2, see idea 2).
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
- `HybridAlgebraicMesh::Theta()` and `CouplingValue()` extract dense blocks with `SparseMatrix::block()`, which scans
  the row from its start: the cost grows with the block size, i.e. with k.
- `BlockJacobi::Setup()` factorizes all the diagonal blocks, while the prolongations 4 and 6 only need those of the
  removed faces (0.02 s at k=0, 0.09 s sequential at k=2).
- Memory: peak 1.33 GB for Cube-tet k=0 n=32 (whole run). The paper's n=64 extrapolates to 10-11 GB, too close to
  13 GB. Running the paper sizes here would first need a memory profile (probably the mesh and the assembly rather
  than the AMG: unverified).
