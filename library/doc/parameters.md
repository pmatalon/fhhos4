[← The fhhos4 library](../../README_library.md)

# Parameters

The public members of `fhhos4::Solver`, read by the setup (`Tolerance` and `MaxIterations` by `Solve()`), with the
library's defaults (`library/include/fhhos4.h`, which details each one):

| Parameter | Default | Values |
|---|---|---|
| `VectorComponents` | `-1`: required if `Dimension` or a degree is inferred, `1` otherwise | `>= 1`: components of the unknown, interleaved in each block (see [Inputs](system.md#inputs)) |
| `Dimension`, `FaceDegree`, `CellDegree` | `-1`: inferred from the blocks (see [Inputs](system.md#inputs)) | `2` or `3`; k `>= 0`; l = k - 1, k or k + 1 |
| `PolynomialSpace` | `PolynomialSpace::Inferred`: from the blocks | Space of the cell and face bases: <ul><li>`P`: degree at most k</li><li>`Q`: degree at most k in each variable (quadrilaterals, hexahedra)</li></ul> |
| `SingularMatrix` | `false` | `true`: the kernel is the near null space given: the constant of each component (pure Neumann), or all the functions given beyond them (the rigid body modes: pure traction). `SetupSym()` only; checked at the setup. The solution is defined up to that kernel |
| `Krylov` | `FCG` if symmetric, `GCR` if non-symmetric | <ul><li>`FCG`: flexible conjugate gradient</li><li>`GCR`: restarted GCR, any matrix</li><li>`None` (`None_` in Python): multigrid cycles alone (e.g. [as a preconditioner](cpp.md#use-as-a-preconditioner))</li></ul> |
| `Restart` | `1` for FCG, `10` for GCR | `>= 1`: search directions stored (two face vectors each). Fewer save memory; GCR may then need more iterations |
| `Tolerance` | `1e-8` | Relative residual of the condensed system at which `Solve()` stops |
| `MaxIterations` | `200` | |
| `Cycle` | `Cycle::K` | <ul><li>`V`</li><li>`W`</li><li>`K`: non-linear: a Krylov method around it must be flexible (FCG, GCR, FGMRES)</li></ul> |
| `PCoarsening` | `PCoarsening::Minus2` | Degrees of the p-levels, from `FaceDegree` k: <ul><li>`Minus2`: k, k-2, ..., 0</li><li>`Minus1`: k, k-1, ..., 0</li><li>`ToZero`: k, 0, for bases not ordered by degree (more iterations from k = 3)</li></ul> |
| `PreSmoothingIterations`, `PostSmoothingIterations` | `2` for Chebyshev (its degree), else `1` | Any value `>= 0` (`0`: none) |
| `PreSmoother` | `Smoother::Default`: `BlockChebyshevFourthKind` if symmetric, else `L1HybridBlockGaussSeidel` (`L1BlockGaussSeidel` with `ReproducibilityMode::OnAnyNThreads`) | One block per face: <ul><li>`BlockGaussSeidel`, `BackwardBlockGaussSeidel`, `SymmetricBlockGaussSeidel`: sequential</li><li>`L1BlockGaussSeidel`, `L1BackwardBlockGaussSeidel`: sequential in each process, l1 Gauss-Seidel on the rows coupled to other processes (see below): convergent whatever the number of processes, the same results on any number of threads</li><li>`HybridBlockGaussSeidel`, `HybridBackwardBlockGaussSeidel`: parallel</li><li>`L1HybridBlockGaussSeidel`, `L1HybridBackwardBlockGaussSeidel`: parallel, l1 Gauss-Seidel (the diagonal of the rows coupled to other threads or processes increased by `L1Theta` times the sum of the \|a_ij\| of these couplings): convergent whatever the numbers of threads and processes</li><li>`DampedBlockJacobi`: damping 2/3</li><li>`BlockChebyshevFirstKind`, `BlockChebyshevFourthKind`: parallel, symmetric matrices only</li></ul> |
| `PostSmoother` | `Smoother::Default`: the adjoint of `PreSmoother` (symmetric cycle) | Same as `PreSmoother` |
| `L1Theta` | `0.5` | `>= 0`: the factor of the l1 terms of the `L1...` smoothers (the sum of the \|a_ij\| of the couplings to other threads, for the hybrid ones, or processes added to the diagonal). From `0.5`, convergent for a symmetric matrix whatever the partition; smaller values may save iterations where the plain sweeps converge (`0`: `BlockGaussSeidel`, `HybridBlockGaussSeidel`) |
| `Prolongation` | `Prolongation::HarmonicExtension` | Values of the faces interior to an aggregate: <ul><li>`HarmonicExtension`: solve their rows of the condensed matrix (shorter setup)</li><li>`Reconstruction`: the article's, from the coarse cell, then one damped Jacobi iteration</li></ul> |
| `KCycleMinimumOperatorReduction` | `2` | `>= 0`. K-cycle: each coarse level runs one iteration of a Krylov method preconditioned by the cycle, and a second one if the first does not reduce the residual by 4. With this parameter, never a second one on a level whose operator is less than this many times smaller, in bytes, than the finer level's (two visits of such a level cost more than the finer level). `0`: no such rule |
| `CoarseningFactor` | `3.8` | Ratio of the numbers of cells of two successive h-levels |
| `CoarseMatrixMaxRows` | `1000` | Maximum number of rows of the coarsest matrix (on several processes: in total) |
| `MinRowsPerProcess` | (MPI only) `50` | On several processes, the coarsening also stops at this many rows per process on average (`0`: only `CoarseMatrixMaxRows`). Each process aggregates its own cells: on many processes, the coarsest matrix cannot get down to `CoarseMatrixMaxRows` |
| `CoarseSolver` | `Cholesky` if symmetric, `LU` if non-symmetric | <ul><li>`Cholesky`, `LU`: Eigen's on one process, MUMPS's on several (in the reproducible modes, Eigen's on the process of rank 0: the coarse matrix gathered there, a serial step)</li><li>`MUMPS`: MUMPS's on one process too; refused on several processes in the reproducible modes</li></ul> |
| `MixedPrecision` | `true` | Levels stored and cycled in single precision (half the memory, faster); the Krylov method, residuals and solution stay in double. Chebyshev and Gauss-Seidel smoothers (all but `DampedBlockJacobi`); otherwise ignored |
| `BlockCSRStorage` | `false` | Level operators stored by blocks of faces when `FaceDegree >= 1` (faster smoothing). Always so in mixed precision; in double, costs a copy of the condensed matrix. Chebyshev and Gauss-Seidel smoothers |
| `Threads` | `0`: OpenMP's setting, at most one per physical core unless `OMP_NUM_THREADS` or OpenMP's binding is set (shared among the processes on the same CPUs) | Number of threads |
| `Reproducibility` | `ReproducibilityMode::None` (`None_` in Python) | The same solution, bit for bit, from run to run, given the same binary, machine, matrix, parameters and number of processes (each process aggregates its own cells); the parallel algorithms that would break it replaced: <ul><li>`None`: no guarantee, the fastest (on several processes, MUMPS's coarse solve varies at the rounding level)</li><li>`OnFixedNThreads`: on the same number of threads per process. On several processes, the coarse matrix gathered on the process of rank 0 (`CoarseSolver::MUMPS` refused)</li><li>`OnAnyNThreads`: whatever the number of threads. Also the default smoother for a non-symmetric matrix `L1BlockGaussSeidel`, the hybrid smoothers refused</li></ul> Not covered: MPI's reductions, deterministic in practice, not by the MPI standard |
| `Comm` | (MPI only) `MPI_COMM_NULL`: one process, no MPI call (MPI may be uninitialized, except with `CoarseSolver::MUMPS`) | <ul><li>`MPI_COMM_SELF`: one process</li><li>the communicator of a distributed system</li></ul> |
| `Verbosity` | `0`: silent | <ul><li>`1`: one line per setup and per solve</li><li>`2`: all messages, with the precision of the levels (and why if `MixedPrecision` is not applied) and the residual at each iteration</li></ul> |

## Constraints

- **Eigen 3.4 or 5.x** in your code, for the header `<fhhos4>`. No Eigen object crosses the boundary of the library:
  its compiled functions take arrays (views), which the header builds from your Eigen matrices and vectors, with your
  version of Eigen and your flags (any version, any alignment settings).
- **Not thread-safe**: one `fhhos4::Solver` object must not be used by several threads at the same time. The library
  parallelizes its own work with OpenMP.
