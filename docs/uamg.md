# U-AMG: how it works in the code

Read before changing `Solver/Multigrid/UncondensedAMG/`, `StaticCondensation.h` or the smoothers. Terms:
`GLOSSARY.md`. Decisions: `docs/adr/` (0005-0010, 0017-0019). Measurements: `PERFORMANCE.md`.

## Static condensation by cells (setup from the blocks)

- Computes the rows of `A` in the pattern of `A_FF`, completed (and sorted) only in the rows where `A_FF` lacks
  couplings of the faces of a cell (zeros dropped by an assembly).
- Keeps the factorizations of the blocks of `A_TT` (ADR-0009) and `A_TT^{-1} A_TF`, and reads `A_FT` at each solve: the
  transpose of `A_TF` (symmetric), or the caller's in place (non-symmetric).
- The products of a row of `A_FT W`: the columns of each cell gathered into a dense vector, its rows of `W` summed
  there, then scattered back (the same sums in the same order, bit for bit). `A_FF` is read from its lower part: only
  its strictly lower part is transposed (the clients give it in full).
- What remains is memory traffic (`A` written as a new matrix, `W` read once per face of each cell): reordering the
  rows, prefetching and register blocking gave nothing.

## Boundary faces

Faces of a single cell are removed by the coarsening (ADR-0019). Program option `-uamg-bfc`
(`UncondensedAMG::CollapseBoundaryFaces`, not in the library): kept, those of an aggregate collapsed into one coarse
face of the aggregate alone (`CoarseFaceAggregates` (k, -1)), numbered after its interfaces, owned by the aggregate's
process. Both in the flat passes and in `HybridAlgebraicMesh::Coarsen()`, bit for bit the same (`CoarseningPassTest`
`*_mixed*`).

Singular systems (`UncondensedAMG::SingularMatrix`, the library's `Solver::SingularMatrix`, set by the program for
`-tc fullneumann`; symmetric): the kernel is the interpolation of 1, z = c on the first DoF of each face
block, 0 on the others (`SingularSolver::Kernel()`). No new input: the hierarchical bases with a constant phi_0 make the
higher-order coefficients of 1 zero. On each level, the kernel is the level's own interpolation of 1 (1 on the
h-levels, where the constant has coordinate 1: P z_c = z_f).
- The coarsest level: `SingularSolver` wraps the coarse solver, its matrix pinned on the DoF of largest |z_i| (row and
  column of the identity) before the factorization; each solve projects its right-hand side onto Im(A), zeroes the
  pinned DoF and projects the solution. Without it, the factorization meets a pivot of the order of the rounding errors,
  and the kernel component of the coarse corrections ruins FCG. The kernel is known: no threshold to detect a null
  pivot (MUMPS's `ICNTL(24)`/`CNTL(3)` would need one; Eigen's `SimplicialLDLT` and `SparseLU` only stop on an exactly
  zero pivot), and the same code for all the coarse solvers. On several processes: the DoF of largest |z_i| on the
  smallest rank that has it; its owner pins its row and its column of `Diag`, the processes that have it as a ghost
  zero its column of `Offd`; the dot products are `Parallel::Dot`.
- The finest level: the output of each cycle projected orthogonally to z (`UncondensedAMG::ExecuteOneIteration()`): the
  smoothers' kernel components would accumulate in the iterates until the rounding errors of A x stall the true
  residual at small tolerances.
- No projection on the intermediate levels: their residuals stay in Im(A) up to rounding (z^T A = 0, R = P^T), so the
  levels in single precision are kept. (The geometric multigrid's per-level hooks `ApplyZeroMeanCondition`/
  `EnforceCompatibilityCondition` are in double only.)
- The library checks |A z| / (|A| |z|) < 1e-8 at the setup (an error otherwise), and projects the condensed right-hand
  side of each solve onto Im(A). The solution is returned up to a multiple of z: the client shifts it (the library has
  no integrals of the basis functions to compute a mean).
- Pinning a DoF in the client instead does not work through the uncondensed interface: it zeroes a column of `A_TF`,
  whose face then has no cell once the p-level extraction drops the zeros.

## p-levels

`P_CS` Minus2 (k, k-2, ..., 0); the operators by block extraction (first coefficients of each block,
`SetupOperatorByBlockExtraction()`), prolongation = injection of the low modes: the bases must be hierarchical.
`PCoarsening::ToZero` (`-p-cs =0`, k then 0) for bases not ordered by degree (more iterations from k = 3).

## h-levels

- Multiple pairwise aggregation (factor 3.8), interface collapsing.
- Interpolation of 1: the constant modes are rescaled before the h-coarsening (`UncondensedLevel::CoarsenMesh()`).
- Coarsening prolongation, program default `ReconstructSmoothedTraceOrInject` (P_F, the papers'): trace of
  `-A_TcTc^{-1} A_TcFc` on the removed faces, smoothed by block Jacobi omega = 2/3; chained as multigrid prolongation.
- Library default `HarmonicExtensionOrInject` (`-coarsening-prolong 10`, `Prolongation::HarmonicExtension`): the
  removed faces of an aggregate by the local solve of their rows of S, the kept faces injected (the paper's
  Q_F^smooth, or P_F with omega = 1; several removed faces: a dense LU, `CoarseningPass::ComputeProlongation()`). No
  cell value: no Theta, no coarse `A_TT` (null on the coarse levels, not rescaled nor extracted for the p-levels:
  `UncondensedLevel::UsesCellBlocks()`); flat passes and chained only (`CheckHarmonicExtension()`, at the setup). The
  cells' interpolation of 1 still rescales `A_TF` for the couplings (the aggregates of P_F), and the faces' stays needed
  (`Q_F`'s injections). The rows of S on the removed faces are added directly into the local system (no merge of the
  rows, any order of the sums).
- Transfers by the prolongations of the passes: on the levels of the flat passes, on one process, whose faces have
  blocks of at least 2 DoFs (`UncondensedLevel::PassProlongationsMinFaceBlockSize`), the prolongations of the passes are
  kept and applied one after the other (and their transposes), without their product P. At degree 0 (all the h-levels
  by default) the factors make the transfers slower (many one-non-zero rows): P. `KeepPassProlongations = false`: P
  everywhere (the reference of `CoarseningPassTest`).

## Non-symmetric matrices

`UncondensedLevel::Symmetric` false, `A_F_T` given (ADR-0007): couplings from the averages of the traces of
`A_TF(T, F)` and `A_FT(F, T)`; reconstruction `-A_TT^{-1} A_TF` by LU of the blocks; smoothing of the prolongations
with the rows of the matrix; Galerkin coarse operators `P^T A P` (the local Galerkin products compute full matrices);
defaults GCR and LU (FCG and Cholesky refused); K-cycle with GCR. Downwind order (`Solver/Multigrid/DownwindOrder.h`;
Bey & Wittum 1997, Hackbusch & Probst 1997): faces sorted after those they are strongly and one-sidedly coupled to,
cycles broken; natural order without one-sided couplings; on a level with one-sided couplings all the Gauss-Seidel
sweeps (backward, symmetric included, pre- and post-smoothing) follow it, Jacobi and SOR ignore it.

## Mixed precision

ADR-0010. `Solver::MixedPrecision` (library default true), `Multigrid::MixedPrecision` (program `-uamg-mixed`, off).
- The setup is the double one (blocks, couplings, local matrices of the operators, prolongations of the passes,
  operators of the p-levels), its results stored in float as produced: the coarse operators assembled in float
  (`LocalMatrices::Assemble<T>()`), the transfers and the p-levels' operators rounded (`SparseMatrixOps::Cast`), the
  smoothers set up from the float operators (`BlockChebyshev`, `BlockDiagonalSolver`: factorizations in double, stored
  in float).
- `UncondensedLevelOf<T>` (T: the precision of the cycle's matrices; the non-template base `UncondensedLevel` holds the
  options, the blocks and `LocalOp`), `Multigrid::CycleIn<T>`; the double instantiation is bitwise the former code.
- The only copy: the finest operator in float (`UncondensedAMG::_singleA`); the p-levels' operators in double are
  released after the setup, the coarsest copied in double for the coarse solver (`Level::SingleOperatorMatrix`;
  `OperatorMatrix` on the finest and coarsest levels only).
- Out of the cycle, in double: `A` (A x for the outer Krylov method: a parallel full product, FCG's own being
  sequential), the residuals and x: each cycle corrects the residual of the current x (iterative refinement:
  `Krylov::None`, non-zero initial guesses).
- Supported (`UncondensedAMG::MixedPrecisionApplies()`): a symmetric matrix on one process, the flat passes with the
  chained prolongations and the Galerkin operators, no option needing `A_F_F` (p-levels after h-levels), the Chebyshev
  smoothers, no zero-mean/compatibility conditions nor exports; otherwise the levels are in double (the reason logged).

## Block CSR storage

`Utils/BlockSparseMatrix.h`; `Solver::BlockCSRStorage` in double (default false), always with the mixed precision;
program `-uamg-bcsr`.
- The operator of the cycle of the levels whose faces have several DoFs (but the coarsest) stored by dense blocks,
  column-major, one column index per block (`Level::BlockOperator`/`SingleBlockOperator`, built by
  `UncondensedLevelOf::OnStartSetup()` from the operator in double; in mixed precision instead of the CSR copy in float,
  in double a copy). Used by the Chebyshev smoothers (`BlockChebyshev::PassBlocks()`, Lanczos) and
  `Level::ResidualOf()/Apply()`.
- Block sizes per block row (`BlockStart`: variable degrees, a future need), all the same today; kernels templated on
  the block size, instances for 2, 3, 4, 5, 6, 9, 10, 15, 16 (`DispatchBlockSize()`, once per pass) and a run-time-size
  one (any, variable).
- The passes sum in the order of the CSR ones; Lanczos' products differ (Eigen 5's sparse product sums with two
  accumulators): results equal to rounding, the same iterations.
- Inverses of the diagonal blocks in the block smoothers (`DiagBlockSolveMethod::Inverse`, computed in double;
  `ApplyInverseBlock()`, fused in the passes): Chebyshev, Gauss-Seidel (hybrid too), SOR, Jacobi
  (`BlockDiagonalSolver::SmootherMethod`); the block Jacobi smoothing of the coarsening's prolongations keeps the LU
  solves (`DefaultMethod`); `-diag-block-solve lu|inv` forces both. Never in the static condensation (ADR-0009).

## Krylov methods

- GCR: restarted every 10 iterations, Notay's economical version (ETNA 37, 2010, Algorithm 3.1); the K-cycle then
  iterates with GCR. BiCGSTAB is not flexible: the setup refuses the K-cycle with it.
- FCG replaces its recurrence residual by the true `b - A x` every 20 iterations (library; 0 in the program) and when it
  meets the tolerance (GCR: at each restart and at the tolerance), then continues if the true one does not: the
  tolerance applies to the true residual. Near the attainable accuracy the recurrence drifts; a tolerance below it runs
  `MaxIterations`, not converged.

## Smoothers

- Default: block Chebyshev of the fourth kind, degree 2, for a symmetric matrix (ADR-0017;
  `Smoother::BlockChebyshevFourthKind`, `bcheb4` in the program; `Pre/PostSmoothingIterations` -1: degree 2), also in
  the program (`SmootherType::Default` resolved by `ApplyProgramArgumentDefaults()`: `bcheb4`, K(2,2) without
  `-cycle`); `bgs,rbgs` for `-uamg-nonsym` and the program's other multigrids.
- Chebyshev (`bcheb4`, `bcheb`; `Solver/FixedPoint/BlockChebyshev.h`): symmetric matrices only (the library refuses
  them otherwise); `lambdaMax` of `D^-1 A` by 10 Lanczos iterations at the setup.
- Non-symmetric matrix on several threads: hybrid block Gauss-Seidel (`hbgs`, `hrbgs`; ADR-0018; scalar kernels at
  k = 0; `Solver/FixedPoint/HybridChunks.h`), sequential sweeps on a level where more than 20% of the non-zeros couple
  rows of different threads (`HybridChunks::MaxExternalFraction`; e.g. GMSH's numbering on the tetrahedral meshes).
