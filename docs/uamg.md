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
block, 0 on the others (`SingularSolver::Kernel()`); with n components, n vectors z_k, c on the DoF k of each block (see
Vector-valued unknowns). No new input: the hierarchical bases with a constant phi_0 make the
higher-order coefficients of 1 zero. With near null space functions beyond the constants (`FaceNearNullSpace`, the
rigid body modes: pure traction elasticity), the kernel is these functions, on the DoFs of degree <= 1 of each face and
0 beyond, orthonormalized (modified Gram-Schmidt, twice: a rotation is not orthogonal to the translations). The kernel
of the finest level is computed once by `UncondensedAMG`'s setup (`Kernel()`), which the library reads for its check
and its projections. On each level, the kernel is the level's own interpolation of 1 (1 on the h-levels, where the
constant has coordinate 1: P z_c = z_f); with the functions, the m constants of their coordinates on their p-level and
the h-levels under it (each function has the coordinate 1 on every block), the functions themselves on the levels in
the given bases.
- The coarsest level: `SingularSolver` wraps the coarse solver, its matrix pinned on one DoF per vector of the kernel
  (row and column of the identity) before the factorization: for each vector in turn, minus its multiples of the
  previous ones that make it zero on their pinned DoFs, the DoF of its largest |coefficient|, so that the kernel
  restricted to the pinned DoFs is triangular with a non-zero diagonal (the pinned matrix is regular); on vectors of
  disjoint supports (the constants), the DoF of largest |z_i| of each. Each solve projects its right-hand
  side onto Im(A), zeroes the pinned DoFs and projects the solution. Without it, the factorization meets a pivot of the order of the rounding errors,
  and the kernel component of the coarse corrections ruins FCG. The kernel is known: no threshold to detect a null
  pivot (MUMPS's `ICNTL(24)`/`CNTL(3)` would need one; Eigen's `SimplicialLDLT` and `SparseLU` only stop on an exactly
  zero pivot), and the same code for all the coarse solvers. On several processes: the DoF of largest |coefficient| on
  the smallest rank that has it, whose coefficients the eliminations take from its owner; its owner pins its row and
  its column of `Diag`, the processes that have it as a ghost zero its column of `Offd`; the dot products are
  `Parallel::Dot`.
- The finest level (`UncondensedAMG::ExecuteOneIteration()`): the right-hand side of each cycle (the residual of the
  outer Krylov method) projected onto Im(A), its output orthogonally to z. The residual regains a kernel component, the
  rounding of A x (about eps |A| |x|), which stays as the residual decreases: sent down the levels, the cycles push
  their iterates along the kernel, and the step lengths of the K-cycle go wrong (a stall at small tolerances on fine
  meshes, pure traction elasticity). Without the projection of the output, the smoothers' kernel components would
  accumulate in the iterates until the rounding errors of A x stall the true residual.
- The other levels in double (`Multigrid::EnforceCompatibilityCondition`, set with `SingularMatrix`): their right-hand
  side, the restricted residual, projected onto Im(A) (`UncondensedLevel::EnforceCompatibilityCondition()`, the kernel
  of the level computed at the setup as the coarsest's: `UncondensedLevel::Kernel`). In exact arithmetic it is in Im(A)
  already (z_c^T P^T r = (P z_c)^T r = 0, P mapping the kernel of the level into that of the finer one); this removes
  the rounding errors of the finer levels. The levels in single precision have none (the hook is in double): their
  right-hand sides come from the finest, projected. Their iterates are not projected (the geometric multigrid's
  `ApplyZeroMeanCondition`): with compatible right-hand sides, nothing pushes them along the kernel.
- The library checks |A z| / (|A| |z|) < 1e-8 for each vector of the kernel after the setup of the multigrid (an error
  otherwise: a Dirichlet face in the system, e.g. a clamped one), and projects the condensed right-hand side of each
  solve onto Im(A). The solution is returned up to an element of the kernel: the client shifts it (the library has no
  integrals of the basis functions to compute a mean).
- Pinning a DoF in the client instead does not work through the uncondensed interface: it zeroes a column of `A_TF`,
  whose face then has no cell once the p-level extraction drops the zeros.

## p-levels

`P_CS` Minus2 (k, k-2, ..., 0); the operators by block extraction (first coefficients of each block,
`SetupOperatorByBlockExtraction()`), prolongation = injection of the low modes: the bases must be hierarchical.
`PCoarsening::ToZero` (`-p-cs =0`, k then 0) for bases not ordered by degree (more iterations from k = 3).

## h-levels

- Multiple pairwise aggregation (factor 3.8), interface collapsing. The aggregation order (breadth-first, from the
  cell that is the strong neighbour of the fewest cells) runs on the threads a level of the search at a time, the same
  order as on one thread (`CoarseningPass::ParallelAggregationOrder()`); the greedy pairing that follows it is
  sequential, and on several threads runs on the cells renumbered in that order (in parallel), whose neighbours are
  then mostly in the cache.
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
  everywhere (the reference of `CoarseningPassTest`). P is chained from the coarse side, P_1 (P_2 (... P_n))
  (`UncondensedLevel::ChainFromCoarseSide()`, on several processes and in the mesh passes too): one product with the rows of the level's
  faces, the others with those of the coarser passes.

## Vector-valued unknowns

`UncondensedAMG::VectorComponents` (`fhhos4::Solver::VectorComponents`, the library only: the program has no
vector-valued discretization) = n: each cell and face block holds the scalar basis times each unit vector, interleaved (phi_i e_c
is the DoF n i + c), so that the first n dim(P^d) DoFs of a block span (P^d)^n: the block sizes of all the levels are
n times the scalar ones (`CreateCoarseLevel()`), and the p-levels extract the first DoFs of each block as in the scalar
case. The constants are the first n DoFs of each block: the rescaling by the interpolation of 1 (one coefficient per
block, that of the scalar basis) applies to them (`UncondensedLevel::VectorComponents`, `ScaleConstantDoFs()`,
`ScaleConstantRows()`). The couplings of the aggregation are the traces of the blocks of `A_TF`, summed over the
components; the h-levels have n DoFs per face (the block smoothers, block CSR). The coarse levels hold the n constants
(the translations), not the rotations, unless they are given (next section). `SingularMatrix`: the kernel of the
n constants, a pinned DoF for each on the coarsest level (`SingularSolver`): the pure Neumann vector Laplacian; with
the rotations given, pure traction elasticity (next section). `Prolongation::Reconstruction`: the n rows of the constants of the
reconstruction of the coarse cell go to the n constants of each removed face (`CoarseningPass::ComputeProlongation()`;
coupled components: coupled rows). Checked by `UAMGLibrary.VectorComponents` and `SingularVectorComponents`: the blocks
⊗ I_n give the scalar iterations and solution on each component, with both prolongations (`CoarseMatrixMaxRows` times
n: it counts rows); `CoupledComponents`, the blocks ⊗ C: the direct solution; on several processes,
`DistributedLibrary.VectorComponents` (the blocks of each process ⊗ I_n, the ghosts unchanged; singular too).

## Near null space beyond the constants (elasticity)

`UncondensedAMG::CellNearNullSpace`, `FaceNearNullSpace` with m columns (the library only, `fhhos4::Solver` with more
near null space functions than components: the rigid body modes, m = 3 in 2D, 6 in 3D): the coefficients of these
functions on the DoFs of degree <= 1 of each cell and face (`UncondensedLevel::FaceNearNullSpaceDoFs` per face), B_T
and B_F. They take the place of the n constants:
- the finest level and the p-levels of degree >= 1 are in the given bases, extracted as usual (k, k - 2, ...; the
  functions are of degree <= 1, within the first DoFs of their blocks);
- the last p-level (degree 0 in the code, `CreateCoarseLevel()`) has m DoFs per block, the coefficients of the
  functions: its prolongation is block-diagonal, B_F on the first DoFs of each face instead of the injection
  (`ProlongIn()`, `RestrictIn()`), and its blocks are the Galerkin ones, B_i^T A_ij B_j
  (`ExtractCoarseMatrix()` with the functions; for `A_TF`, B_T^T A_TF B_F). There each function has the coordinate 1
  on every cell and face, as the constant of a scalar unknown after the rescaling by its interpolation of 1. On
  several processes, each one gives the functions of its own cells and faces; those of its ghost faces, which the
  columns of `Offd` and of `A_TF` read, come from their owners (`GhostFunctions()`, a halo exchange per function, in
  `SetupOperatorByBlockExtractionDistributed()`);
- the h-levels have m DoFs per face, and the passes treat them as constants (`UncondensedLevel::VectorComponents` and
  `CoarseningPass::VectorComponents` = m): the injection of the kept faces copies the same m coordinates, i.e. the
  same global function, to all the fine faces of a coarse face.

Precision: with the mixed precision, the finest level and the p-levels of degree >= 1 are in single precision, the
level of the functions and the h-levels under it in double (`CreateCoarseLevel()`): there a rotation is a global
function, whose energy is the small difference of large coefficients (those of its translation part, as |x|^2 / h^2),
lost in single precision (the cycle diverges on a fine mesh). The cycle changes precision at that step
(`Multigrid::CycleIn()`: the restricted residual cast to double, `CoarseCorrection<double>()`, its result cast back),
as it does for the coarsest solve; the level in double extracts its blocks from the operators in double of the finer
level (`SetupOperatorByBlockExtraction()`).

`SingularMatrix` with the functions (pure traction): the kernel is the functions (see Singular systems, under
Boundary faces); on the level of the functions and the h-levels under it, the m constants, which the passes preserve
(P z_c = z_f, as for the constants of the components), so that the coarsest level, pinned on m DoFs, has the kernel
of the m constants. The cells of the functions must be those of the kernel (-A_TT^{-1} A_TF on their faces, what the
interpolation of a rigid motion gives in HHO): otherwise the Galerkin blocks of the level of the functions lose the
kernel.

Keeping the degree 1 on the h-levels without them (`-hp-cs h`) does not carry the rotations: the injection gives all
the fine faces of a coarse face the same coefficients, and the constant part of a rotation differs from one face to
another. Checked by `UAMGLibrary.NearNullSpaceBeyondTheConstants` (algebraic: a third function on 2 components) and
`SingularNearNullSpace` (a kernel of 3 functions, not orthogonal), on several processes by
`DistributedLibrary.NearNullSpaceBeyondTheConstants` (also singular); elasticity itself only on Disk++'s
`linear_elasticity` (`docs/clients.md`, one process: Dirichlet, a clamped beam, pure traction).

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
  smoothers set up from the float operators (`BlockChebyshev`, `BlockGaussSeidel`, `GaussSeidel`,
  `BlockDiagonalSolver`: factorizations in double, stored in float; the l1 terms, and the downwind order of the
  h-levels, from the float values).
- `UncondensedLevelOf<T>` (T: the precision of the cycle's matrices; the non-template base `UncondensedLevel` holds the
  options, the blocks and `LocalOp`), `Multigrid::CycleIn<T>`; the double instantiation is bitwise the former code.
- The only copy: the finest operator in float (`UncondensedAMG::_singleA`); the p-levels' operators in double are
  released after the setup, the coarsest copied in double for the coarse solver (`Level::SingleOperatorMatrix`;
  `OperatorMatrix` on the finest and coarsest levels only).
- Out of the cycle, in double: `A` (A x for the outer Krylov method: a parallel full product, FCG's own being
  sequential), the residuals and x: each cycle corrects the residual of the current x (iterative refinement:
  `Krylov::None`, non-zero initial guesses).
- The smoothers in single precision (`IterativeSolver::SolveInSinglePrecision()`, the storages `_singleMatrix`...):
  the vectors and matrices in float, the sums in double, rounded when stored. Gauss-Seidel (`BlockGaussSeidel` and its
  hybrid and l1 variants; `GaussSeidel` on the levels of degree 0): its kernels templated on the precision; the block
  solves by `BlockDiagonalSolver::ApplyInverse()` on the right-hand side in double (the LU solves in float with
  `-diag-block-solve lu`).
- Supported (`UncondensedAMG::MixedPrecisionApplies()`): symmetric and non-symmetric matrices (on several processes
  too: `docs/mpi.md`), the flat passes with the chained prolongations and the Galerkin operators, no option needing
  `A_F_F` (p-levels after h-levels), the Chebyshev and Gauss-Seidel smoothers (of the block size of the finest level
  and scalar), no zero-mean condition nor exports (the compatibility condition of `SingularMatrix` runs on the levels
  in double only); otherwise the levels are in double (the reason logged).

## Block CSR storage

`Utils/BlockSparseMatrix.h`; `Solver::BlockCSRStorage` in double (default false), always with the mixed precision;
program `-uamg-bcsr`.
- The operator of the cycle of the levels whose faces have several DoFs (but the coarsest) stored by dense blocks,
  column-major, one column index per block (`Level::BlockOperator`/`SingleBlockOperator`, built by
  `UncondensedLevelOf::OnStartSetup()` from the operator in double; in mixed precision instead of the CSR copy in float,
  in double a copy). Used by the Chebyshev smoothers (`BlockChebyshev::PassBlocks()`, Lanczos), the block Gauss-Seidel
  (`BlockGaussSeidel::BlockRowRhs()`; the hybrid chunks by ranges of blocks, `HybridChunks::Setup(BlockSparseMatrix)`)
  and `Level::ResidualOf()/Apply()`. On several processes, `Diag` (`Offd`, whose columns are ghosts, stays in CSR).
- Block sizes per block row (`BlockStart`: variable degrees, a future need), all the same today; kernels templated on
  the block size, instances for 2, 3, 4, 5, 6, 9, 10, 15, 16 (`DispatchBlockSize()`, once per pass) and a run-time-size
  one (any, variable).
- The passes and sweeps sum in the order of the CSR ones (the Gauss-Seidel sweeps: the results of CSR, bit for bit);
  Lanczos' and the residuals' products differ (Eigen 5's sparse product sums with two accumulators): results equal to
  rounding, the same iterations.
- Inverses of the diagonal blocks in the block smoothers (`DiagBlockSolveMethod::Inverse`, computed in double;
  `BlockDiagonalSolver::ApplyInverse()` in Chebyshev, fused in the passes, and in Gauss-Seidel in single precision):
  Chebyshev, Gauss-Seidel (hybrid too), SOR, Jacobi
  (`BlockDiagonalSolver::SmootherMethod`); the block Jacobi smoothing of the coarsening's prolongations keeps the LU
  solves (`DefaultMethod`); `-diag-block-solve lu|inv` forces both. Never in the static condensation (ADR-0009).

## Krylov methods

- GCR: restarted every 10 iterations, Notay's economical version (ETNA 37, 2010, Algorithm 3.1); the K-cycle then
  iterates with GCR. FCG keeps 1 direction (FCG(1), truncated). The library's `Restart` sets both (the outer method
  only, not the K-cycle's).
- K-cycle (`Multigrid::KrylovForKCycle()`): on each coarse level, one iteration of FCG or GCR preconditioned by the
  cycle, and a second one if the first does not reduce the residual by 4 (Notay). `KCycleMinimumOperatorReduction` R
  (library, default 2; program `-kcycle-min-reduction`, default 0: the papers) adds a criterion known at the setup
  (`Level::OneKCycleIteration`, listed in the setup log): never a second iteration on a level whose operator is less
  than R times smaller than the finer level's, in bytes (`Level::OperatorBytes()`: the values in the level's
  precision, not the indices, which depend on the storage), since the cycle is bound by the memory. Such a level skips the
  residual update and its two norms. One rule for all the systems: it applies where the operators decrease slowly,
  i.e. at the first aggregation in 3D on simplices, and on the level of the near null space functions (in double
  under levels in single precision).
- FCG replaces its recurrence residual by the true `b - A x` every 20 iterations (library; 0 in the program) and when it
  meets the tolerance (GCR: at each restart and at the tolerance), then continues if the true one does not: the
  tolerance applies to the true residual. Near the attainable accuracy the recurrence drifts; a tolerance below it runs
  `MaxIterations`, not converged.
- Vector operations (FCG, GCR, the K-cycle, the Lanczos estimate of the Chebyshev smoothers): by blocks of 4096 rows
  in parallel (`Parallel::RowBlocks`, more rows per block above 4M rows), each block's dot products summed with the
  others' in the order of the blocks (`Parallel::Dot()` too): the results depend on the number of rows, not of threads;
  on one block, Eigen's dot product of the whole vectors, bit for bit. The operations of an iteration are fused into
  passes, each block updated then reduced while in the cache (`Parallel::SumsOverRowBlocks()`): FCG three passes (the
  z.Ad_k; d, A d, d.Ad and r.d; x, r and |r|^2), GCR j + 3 (one per Gram-Schmidt step, each subtracting the previous
  projection and computing the next dot product; c_j normalized with c_j.r; r and |r|^2), the K-cycle one for the
  dot products of an iteration, one for r and |r|^2. Each coefficient goes through the same operations as one
  operation after the other: the fused passes give the results of the separate ones bit for bit.

## Vectors of the solve

The solve allocates no vector: a vector of more than 32 MB (glibc's largest reused block) would be fresh memory at each
allocation, its page faults paid together by the threads of a process (the processes of a run each pay their own).
Each vector the solve uses is allocated at the setup by the component that uses it, and written in place:
- the levels: `Level::CycleVectors` (the residual, its restriction, the coarse correction, the K-cycle's vectors),
  allocated at the end of `Multigrid::SetupLevels()` (`Level::AllocateCycleVectors()`); the transfers write in them
  (`RestrictInto()`, `ResidualInto()`, `ApplyInto()`), and `ProlongAndAdd()` adds P e to x directly (no vector for P e:
  on several processes, Diag e then Offd e_ghosts added to x);
- the smoothers: their work vectors at their setup (`BlockChebyshev::AllocateWorkVectors()`, `GaussSeidel`'s products
  of the sweep, shared with the post-smoother that copies the setup); they write the residual and A x of their last
  iteration in the caller's vectors (`IterativeSolver::SolveInto()`; `Smoother::SmoothAndComputeResidualOrAx()`);
- the multigrid as the preconditioner of the Krylov method writes its A e in the caller's vector (`SolveInto()`,
  `SolverPreconditioner::ApplyAndComputeAe()`): FCG's new direction (its vectors, those of its `Truncation` + 1
  directions, and GCR's `Restart` vectors z_j and c_j, at their setup).

The setup's large arrays (those of the flat passes, `CoarseningPass`, `IndexLists`, `LocalMatrices`, and the
condensation's factors) are `UninitializedVector`s (`Parallelism.h`): `resize()` writes nothing, and a parallel loop
writes them first (`ParallelAssign()` for those that need a value), not `vector<T>(n)`'s zero-fill by the calling
thread. On a large problem each is fresh memory, whose page faults cost ~4.5 times more on one thread than on 8
(`PERFORMANCE.md`). The lists of matrices that grow by `push_back()` (the prolongations of the passes) are deques:
Eigen's `SparseMatrix` move constructor is not `noexcept`, so a vector would copy them at each reallocation.

Exceptions: the coarse solver's solution (the coarsest level: small); the stand-alone solvers (`-s gs`, `-s uamg`...),
whose stopping criteria read the residual of each iteration, keep `IterativeSolver::Solve()`'s loop (a vector per
iteration). The first iteration writes the vectors first: on a NUMA machine their pages land where the solve's threads
write them, not on the node of the setup's thread.

## Smoothers

- Default: block Chebyshev of the fourth kind, degree 2, for a symmetric matrix (ADR-0017;
  `Smoother::BlockChebyshevFourthKind`, `bcheb4` in the program; `Pre/PostSmoothingIterations` -1: degree 2), also in
  the program (`SmootherType::Default` resolved by `ApplyProgramArgumentDefaults()`: `bcheb4`, K(2,2) without
  `-cycle`). For a non-symmetric matrix, l1 hybrid block Gauss-Seidel (`L1HybridBlockGaussSeidel`, `l1hbgs,l1hrbgs`,
  ADR-0018), also for `-uamg-nonsym`; `bgs,rbgs` for the program's other multigrids.
- Chebyshev (`bcheb4`, `bcheb`; `Solver/FixedPoint/BlockChebyshev.h`): symmetric matrices only (the library refuses
  them otherwise); `lambdaMax` of `D^-1 A` by 10 Lanczos iterations at the setup.
- Hybrid block Gauss-Seidel (`hbgs`, `hrbgs`; ADR-0018; scalar kernels at k = 0; `Solver/FixedPoint/HybridChunks.h`):
  sequential sweeps on a level where more than 20% of the non-zeros couple rows of different threads
  (`HybridChunks::MaxExternalFraction`; e.g. GMSH's numbering on the tetrahedral meshes).
- Jacobi between the parts (the chunks of the hybrid sweeps, the processes) is not a contraction at k = 2 in 3D under
  that threshold: the cycle diverges; nor at k = 0 on several processes (Cube-tet n=32, 2 processes: GCR stagnates at
  3.7e-2). l1 Gauss-Seidel (`l1hbgs`, `l1hrbgs`: `BlockGaussSeidel` with `l1`, and `GaussSeidel` with `l1` at k = 0;
  ADR-0018): the rows coupled to other parts get their l1 term, theta sum_j |a_ij| over those couplings
  (`IterativeSolver::L1Terms()`: the rows of `Offd`, `HybridChunks::AddExternalAbsoluteSums()`), on their diagonal (in
  the factorized blocks, refactorized when the chunks change: `SetupL1Terms()`; in the scalar kernels,
  `GaussSeidel::Diagonal()`), and times the old x_i in the right-hand side (`SweepRhs()`). theta: each smoother's
  (`_l1Theta`, 0: no l1 terms), from `Multigrid::L1Theta` through `SmootherFactory::Create()` (library `L1Theta`, program
  `-gs-l1`), 1/2 by default: M + M^T - A >= D for a symmetric A, whatever the parts. No row has an l1 term on one thread
  and one process.
