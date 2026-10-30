// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Pierre Matalon
#pragma once
#include <cstdint>
#include <memory>
#include <type_traits>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "fhhos4_error.h"
#include "fhhos4_version.h"
#if FHHOS4_MPI
#include <mpi.h>
#endif

// fhhos4: the algebraic multigrid of fhhos4 for hybrid discretizations (HHO, HDG...) of scalar
// diffusion problems, the Uncondensed AMG (U-AMG) of
//   D. A. Di Pietro, F. Hülsemann, P. Matalon, P. Mycek, U. Rüde, "Algebraic multigrid preconditioner for statically
//   condensed systems arising from lowest-order hybrid discretizations", SIAM J. Sci. Comput. (2023),
// extended to k >= 1.
// Solves a hybrid system (cell and face unknowns), symmetric positive definite (SetupSymmetric()) or not
// (SetupNonsymmetric(), e.g. convection-diffusion),
//     [ A_TT  A_TF ] [x_T]   [b_T]
//     [ A_FT  A_FF ] [x_F] = [b_F],
// by static condensation: the AMG solves the condensed system A x_F = b, A = A_FF - A_FT A_TT^{-1} A_TF and
// b = b_F - A_FT A_TT^{-1} b_T: SetupSymmetric() or SetupNonsymmetric() with the blocks, then Solve(b_T, b_F, x_T, x_F).
// The library condenses the system and computes the cell unknowns: no static condensation in the calling code.
// See README_library.md.
//
// The Krylov method preconditioned by the AMG is part of the library (Solver::Krylov), which does not give the multigrid
// cycle alone: with the default K-cycle, the AMG is a non-linear preconditioner (it solves the coarse levels with inner
// FCG or GCR iterations), which requires a flexible Krylov method: the flexible conjugate gradient (FCG, the default for
// a symmetric matrix) or GCR.
//
// The inputs are Eigen objects, read in place, never copied nor converted: the sparse matrices in the format of the
// library, fhhos4::SparseMatrix (row-major, with the index type SparseIndex) or an Eigen::Map of one, the vectors as
// Eigen vectors. Another Eigen sparse matrix (column-major, another index type) converts by Eigen, a copy in the
// calling code: fhhos4::SparseMatrix B = A; (a symmetric column-major one: TransposedView(), no copy). Assembling the
// matrices in the format of the library avoids the copies. The helpers at the end of this file extract the blocks
// from the global matrix of the system (copies).
//
// The functions compiled in the shared library fhhos4 take arrays only: the Eigen objects of this header are those of
// the calling code, compiled with its own version of Eigen (the inline functions of this header give their arrays to
// the library). The library never allocates memory that the calling code frees, nor the other way around.
//
// Compiled with MPI (FHHOS4_MPI, -DENABLE_MPI=ON), a Solver works on the processes of its communicator (Solver::Comm;
// without it, on the calling process alone, as without MPI), on a system distributed by rows: each process gives the
// rows of its own cells and faces, the columns of the faces numbered locally, its own faces first, then the other faces
// of its cells (its ghosts), whose global numbers it gives too (GhostNumbers); the same for the columns of the cells of
// A_FT (non-symmetric). See the setups below, and README_library.md.
namespace fhhos4
{
	// The index type of the sparse matrices of the library: int (32 bits), or 64 bits if fhhos4 is compiled with
	// -DSMALL_INDEX=OFF
	using SparseIndex = FHHOS4_SPARSE_INDEX;

	// The sparse matrices of the library: Eigen's, row-major, with the index type SparseIndex, compressed. The setups
	// read them in place, or an Eigen::Map of one (SparseMatrixMap), e.g. of the arrays of another code.
	using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, SparseIndex>;

	// A sparse matrix in the format of the library on arrays it doesn't own (no copy): e.g. TransposedView(),
	// Solver::CondensedMatrix()
	using SparseMatrixMap = Eigen::Map<const SparseMatrix>;

	namespace detail
	{
		// The vectors given to the compiled functions: the arrays of the calling code's vectors, read (or written) in place
		struct ConstVectorView
		{
			const double* Data;
			std::int64_t Size;
		};

		struct VectorView
		{
			double* Data;
			std::int64_t Size;
		};
	}

	// The implementation: the arrays given to the compiled functions, and what they return (the Python module uses them
	// too, on the arrays of NumPy)
	namespace detail
	{
		// A sparse matrix in the CSR format (compressed sparse row), not owned: arrays read in place. The non-zeros of the
		// row i are (ColumnIndices[p], Values[p]) for RowPointers[i] <= p < RowPointers[i + 1], with RowPointers[0] = 0 and
		// the columns strictly increasing in each row (no duplicates). RowPointer and Column: std::int32_t or
		// std::int64_t, independently.
		template <typename RowPointer, typename Column>
		struct BasicCSRView
		{
			static_assert((std::is_same<RowPointer, std::int32_t>::value || std::is_same<RowPointer, std::int64_t>::value) &&
			              (std::is_same<Column, std::int32_t>::value || std::is_same<Column, std::int64_t>::value),
			              "fhhos4: the indices are std::int32_t or std::int64_t");

			std::int64_t Rows = 0;
			std::int64_t Cols = 0;
			const RowPointer* RowPointers = nullptr; // Rows + 1
			const Column* ColumnIndices = nullptr;   // RowPointers[Rows]
			const double* Values = nullptr;          // RowPointers[Rows]

			BasicCSRView() = default;

			BasicCSRView(std::int64_t rows, std::int64_t cols, const RowPointer* rowPointers, const Column* columnIndices, const double* values)
				: Rows(rows), Cols(cols), RowPointers(rowPointers), ColumnIndices(columnIndices), Values(values)
			{}

			std::int64_t NonZeros() const
			{
				return RowPointers ? (std::int64_t)RowPointers[Rows] : 0;
			}
		};

		// The matrices of the setups: the index type of the library
		using CSRView = BasicCSRView<SparseIndex, SparseIndex>;

		// The arrays of an Eigen sparse matrix (compressed, row-major, of doubles, with the indices Index), read in place
		template <typename Index, typename Derived>
		BasicCSRView<Index, Index> View(const Eigen::SparseCompressedBase<Derived>& M)
		{
			if (!M.isCompressed())
				throw Error("fhhos4: the Eigen matrix is not compressed: call its makeCompressed() first.");
			return BasicCSRView<Index, Index>(M.rows(), M.cols(), M.outerIndexPtr(), M.innerIndexPtr(), M.valuePtr());
		}

		// A sparse matrix in the CSR format of the library, owned by the library (allocated and freed by it): the result of
		// the compiled helpers that copy
		class FHHOS4_API CSRMatrix
		{
		public:
			struct Impl; // defined in the library

			CSRMatrix();
			explicit CSRMatrix(std::unique_ptr<Impl> impl);
			~CSRMatrix();
			CSRMatrix(CSRMatrix&&) noexcept;
			CSRMatrix& operator=(CSRMatrix&&) noexcept;
			CSRMatrix(const CSRMatrix&) = delete;
			CSRMatrix& operator=(const CSRMatrix&) = delete;

			CSRView View() const&; // 0 x 0 if empty
			CSRView View() const&& = delete;

		private:
			std::unique_ptr<Impl> _impl;
		};

		// A copy, owned by the calling code: an Eigen matrix
		inline SparseMatrix ToEigen(const CSRMatrix& M)
		{
			CSRView v = M.View();
			return SparseMatrix(SparseMatrixMap(v.Rows, v.Cols, v.NonZeros(), v.RowPointers, v.ColumnIndices, v.Values));
		}
	}

	// On several processes (MPI): the global numbers of the ghosts of a process (the faces, or cells, of its rows owned by
	// other processes), in the order of their local numbers (after those of its own faces, or cells), read in place. The
	// global numbers: the faces of the process of rank 0 first, then those of the process of rank 1, etc., each process's
	// in the order of its local numbers (so that the faces of each process are numbered contiguously). Any order of the
	// ghosts, no repetition.
	struct GhostNumbers
	{
		std::int64_t Count = 0;
		const std::int64_t* Numbers = nullptr; // Count

		GhostNumbers() = default;
		GhostNumbers(std::int64_t count, const std::int64_t* numbers)
			: Count(count), Numbers(numbers)
		{}
	};

	// The Krylov method preconditioned by one cycle of the AMG (Solver::Krylov)
	enum class Krylov
	{
		Default,  // FCG for a symmetric matrix (SetupSymmetric...()), GCR otherwise (SetupNonsymmetric...())
		FCG,      // flexible conjugate gradient (fhhos4 -s fcguamg), for a symmetric positive definite matrix
		GCR,      // GCR(10), generalized conjugate residual restarted every 10 iterations, Notay's economical version
		          // (fhhos4 -s gcruamg). Minimizes the residual norm: for any matrix
		None      // the cycles of the AMG alone, as a stationary iterative method (fhhos4 -s uamg): slower at k >= 1
	};

	// The multigrid cycle (Solver::Cycle)
	enum class Cycle
	{
		V,
		W,
		K // K-cycle: each coarse level is solved by 2 iterations of FCG, or of GCR with Krylov::GCR or a non-symmetric
		  // matrix. Non-linear: it requires a flexible Krylov method (FCG or GCR)
	};

	// The polynomial space of the bases of the cells and faces (Solver::PolynomialSpace)
	enum class PolynomialSpace
	{
		P, // P^k: the polynomials of degree at most k; (k + n)! / (k! n!) functions in n variables
		Q  // Q^k: the polynomials of degree at most k in each variable (tensor products, on quadrilaterals and hexahedra);
		   // (k + 1)^n functions in n variables
	};

	// The p-coarsening (Solver::PCoarsening): the degrees of the p-levels, from FaceDegree down to 0, before the h-levels
	// (aggregation) of degree 0. A p-level of degree d keeps the first functions of each block of the finer level: those
	// of the polynomials of degree <= d (P^d, or Q^d: see PolynomialSpace) if the bases are ordered by degree.
	enum class PCoarsening
	{
		Minus2,    // k, k - 2, ..., 0 (fhhos4 -p-cs -2)
		Minus1,    // k, k - 1, ..., 0 (-p-cs -1)
		DivideBy2, // k, k / 2, ..., 0 (-p-cs /2)
		ToZero     // k, 0 (-p-cs =0): the only p-level keeps the first function of each block, the constant one: for bases
		           // not ordered by degree. Usually more iterations than the others at k >= 3
	};

	// The prolongation of the h-coarsening (Solver::Prolongation): the faces of the boundary of an aggregate of cells take
	// the value of their coarse face; the faces interior to it (removed on the coarse level) get theirs from them:
	enum class Prolongation
	{
		HarmonicExtension, // the default: the solution of their rows of the condensed matrix, the faces of the boundary
		                   // given. No value of the cells, whose bases then only enter the choice of the aggregates; a
		                   // shorter setup (fhhos4 -coarsening-prolong 10)
		Reconstruction     // the paper's: the value of the coarse cell, reconstructed from its faces (the inverse of the
		                   // static condensation), then one damped Jacobi iteration of the condensed matrix on these
		                   // faces (fhhos4 -coarsening-prolong 6, the program's default)
	};

	// The smoothers (Solver::PreSmoother and PostSmoother; fhhos4 -smoothers). The blocks: one per face. The Gauss-Seidel
	// smoothers are sequential, except the hybrid ones. On several processes (MPI), they are Gauss-Seidel in each
	// process (each thread for the hybrid ones), Jacobi between the processes.
	enum class Smoother
	{
		Default,                        // chosen by the setup (see Solver::PreSmoother and PostSmoother)
		BlockGaussSeidel,               // bgs
		BackwardBlockGaussSeidel,       // rbgs
		SymmetricBlockGaussSeidel,      // sbgs: forward then backward
		HybridBlockGaussSeidel,         // hbgs: Gauss-Seidel in the rows of each thread, Jacobi between the threads:
		                                // parallel, except on the levels where the rows of the threads are poorly
		                                // separated; the results depend on the number of threads, not on the run. On 1
		                                // thread: BlockGaussSeidel
		HybridBackwardBlockGaussSeidel, // hrbgs
		DampedBlockJacobi,              // bj23: block Jacobi damped by 2/3. Parallel
		BlockChebyshevFirstKind,        // bcheb: block Chebyshev, a polynomial in D^-1 A (D: the block diagonal) whose
		                                // degree is the number of smoothing iterations (each one, a block Jacobi
		                                // iteration with its own damping factor), minimal on [lambdaMax/10, lambdaMax].
		                                // Parallel, and its results don't depend on the number of threads. Symmetric
		                                // positive definite matrices only
		BlockChebyshevFourthKind        // bcheb4: the same, of the fourth kind (Lottes 2023): no lower bound to choose
	};

	// The solver of the coarsest level (Solver::CoarseSolver). On several processes (MPI), the factorizations are MUMPS's,
	// distributed over the processes.
	enum class CoarseSolver
	{
		Default,  // Cholesky for a symmetric matrix (SetupSymmetric...()), LU otherwise (SetupNonsymmetric...())
		Cholesky, // sparse Cholesky factorization, for a symmetric positive definite matrix
		LU,       // sparse LU factorization
		MUMPS     // MUMPS's sparse factorization: LDL^T for a symmetric matrix (positive definite), LU otherwise. Needs
		          // the library compiled with MUMPS (ENABLE_MUMPS; always with MPI)
	};

	struct Result
	{
		int Iterations = 0;
		double RelativeResidual = 0; // |b - A x| / |b|, in the condensed system
		bool Converged = false;      // RelativeResidual < Tolerance
	};

	// Set the parameters, then call a setup function, then Solve(). The parameters are read by the setup (Tolerance and
	// MaxIterations by Solve()). The defaults are those of the program fhhos4 (-s fcguamg), used in the papers, except the
	// smoothers (PreSmoother, and the number of its iterations), the prolongation (Prolongation) and the precision of the
	// levels (MixedPrecision).
	// Not thread-safe: one object must not be used by several threads at the same time.
	class FHHOS4_API Solver
	{
	public:
		//------------- Discretization: to set -------------//

		int Dimension = -1;  // 2 or 3
		int FaceDegree = -1; // k: degree of the polynomials on the faces
		int CellDegree = -1; // l: degree of the polynomials in the cells (k - 1, k or k + 1 in HHO)
		// The polynomial space of the cells and faces, P^l and P^k, or Q^l and Q^k: their numbers of DoFs, and those of the
		// p-levels
		fhhos4::PolynomialSpace PolynomialSpace = fhhos4::PolynomialSpace::P;
		// The matrix is singular, its kernel spanned by the interpolation of 1 (cellInterpOfOne, faceInterpOfOne: the
		// constant functions), e.g. the pure Neumann problem: SetupSymmetric() only. The setup checks that
		// the interpolation of 1 is in the kernel; each solve projects the condensed right-hand side onto the image of the
		// matrix (the data of a pure Neumann problem are compatible up to the rounding or quadrature errors; at Verbosity
		// >= 1, a message if the component removed is larger). The solution is defined up to a multiple of the
		// interpolation of 1 (a constant): shift it to the normalization you want (e.g. zero mean).
		bool SingularMatrix = false;

		//--------------------- Solver ---------------------//

		// The Krylov method preconditioned by one cycle of the AMG. Default: FCG for a symmetric matrix, GCR otherwise.
		fhhos4::Krylov Krylov = fhhos4::Krylov::Default;
		double Tolerance = 1e-8; // stopping criterion of Solve(): |b - A x| / |b| < Tolerance
		int MaxIterations = 200; // Solve() stops after MaxIterations iterations, converged or not

		//-------------------- Multigrid -------------------//

		fhhos4::Cycle Cycle = fhhos4::Cycle::K;
		fhhos4::PCoarsening PCoarsening = fhhos4::PCoarsening::Minus2; // the degrees of the p-levels (k >= 1)
		// Number of iterations of the pre- and post-smoothers (>= 0, 0: none); -1 (default): the smoother's, 2 for the
		// Chebyshev smoothers (their degree), 1 for the others
		int PreSmoothingIterations = -1;
		int PostSmoothingIterations = -1;
		// Smoothers. PreSmoother Default: chosen by the setup:
		// - for a symmetric matrix (SetupSymmetric...()), BlockChebyshevFourthKind: parallel, and its results don't depend
		//   on the number of threads;
		// - for a non-symmetric matrix (the Chebyshev smoothers require a symmetric one), HybridBlockGaussSeidel on
		//   several threads (parallel), BlockGaussSeidel on 1.
		// PostSmoother Default: the adjoint of PreSmoother, which makes the cycle symmetric: the sweep in the other
		// direction for Gauss-Seidel (BackwardBlockGaussSeidel after BlockGaussSeidel, HybridBackwardBlockGaussSeidel after
		// HybridBlockGaussSeidel, and the other way around), PreSmoother itself if it is symmetric (Chebyshev, Jacobi,
		// SymmetricBlockGaussSeidel).
		// Non-symmetric matrix (SetupNonsymmetric...()): on a level with one-sided couplings (convection), the Gauss-Seidel
		// smoothers keep their type, block size and number of sweeps, but not their direction: all their sweeps
		// (backward and symmetric ones included) follow the downwind order of the level. Jacobi ignores it.
		fhhos4::Smoother PreSmoother = fhhos4::Smoother::Default;
		fhhos4::Smoother PostSmoother = fhhos4::Smoother::Default;
		// The prolongation of the h-coarsening: from the faces of the boundary of each aggregate to its interior faces
		fhhos4::Prolongation Prolongation = fhhos4::Prolongation::HarmonicExtension;
		double CoarseningFactor = 3.8;     // target ratio of the numbers of cells of two successive h-levels
		int CoarseMatrixMaxSize = 1000;    // the coarsening stops when the matrix has at most this many rows
		// Solver of the coarsest level. Default: Cholesky for a symmetric matrix, LU otherwise.
		fhhos4::CoarseSolver CoarseSolver = fhhos4::CoarseSolver::Default;
		// Mixed precision: the matrices of the levels of the AMG (operators, transfers, smoothers) stored and its cycles
		// run in single precision (float), but the coarsest solve: half the memory of the levels, faster cycles. The setup
		// computes in double (the same levels, rounded), and the Krylov method, the matrix A, the residuals and the
		// solution stay in double: the tolerance can be that of the double precision, usually in the iterations of the
		// double precision. Only for a symmetric matrix, on one process, with the Chebyshev smoothers (the defaults);
		// otherwise the levels are in double, as with false.
		bool MixedPrecision = true;
		// The operator of the levels whose faces have several DoFs (FaceDegree >= 1) stored by blocks (block CSR: one
		// column index per block of coupled faces instead of one per coefficient) for the smoothers and the cycle: their
		// passes over the matrix are 1.4-2.5x faster. Always so in mixed precision (MixedPrecision; no copy: instead of
		// the copy in single precision). In double precision (MixedPrecision false), with true: at the cost of a copy of
		// the condensed matrix A (about 2/3 of its memory). The Chebyshev smoothers (the defaults), on one process;
		// otherwise in CSR.
		bool BlockCSRStorage = false;

		//-------------------- Execution -------------------//

		int Threads = 0;   // number of OpenMP threads during the calls of the library; 0: the current OpenMP setting, but
		                   // at most one thread per physical core (no hyper-threading) unless OMP_NUM_THREADS or OpenMP's
		                   // binding (OMP_PROC_BIND, OMP_PLACES) is set. On several processes (MPI): the cores divided
		                   // among the processes of Comm that run on the same CPUs (not bound to distinct cores)
		int Verbosity = 0; // 0: nothing; 1: one line per setup (number of levels, time) and per solve (iterations,
		                   // residual, time); 2: all the messages of fhhos4, with the residual at each iteration.
		                   // On several processes: the messages of the process of rank 0.

#if FHHOS4_MPI
		// The processes of the solver (the library is compiled with MPI). MPI_COMM_NULL (the default: not set) or
		// MPI_COMM_SELF: the system of the calling process, solved by it alone, as without MPI; not set, the solver makes no
		// MPI call (MPI may not be initialized), but CoarseSolver::MUMPS needs MPI initialized; forgotten for a distributed
		// system, the setup refuses its inputs (ghosts, A_FF with more columns than rows) and says to set Comm. A
		// communicator of several processes: one system distributed by rows over them (see SetupSymmetric()); every setup
		// and solve is then collective, called by all the processes of the communicator, with the same parameters.
		MPI_Comm Comm = MPI_COMM_NULL;
#endif

		Solver();
		~Solver();
		Solver(Solver&&) noexcept;
		Solver& operator=(Solver&&) noexcept;
		Solver(const Solver&) = delete;
		Solver& operator=(const Solver&) = delete;

		// The inputs of the setups:
		// - the matrices: Eigen sparse matrices in the format of the library, fhhos4::SparseMatrix (compressed, row-major,
		//   with the index type SparseIndex) or an Eigen::Map of one (SparseMatrixMap), read in place, never copied. Other
		//   Eigen matrices convert by Eigen, a copy: fhhos4::SparseMatrix B = A; (a symmetric column-major one:
		//   TransposedView(A), no copy). The setups check the sizes and structure of the inputs against Dimension and the
		//   degrees, and their CSR arrays (one pass on the indices);
		// - cellInterpOfOne, faceInterpOfOne: the interpolation of the function 1 on the polynomial bases of the cells and
		//   faces, one coefficient per cell and per face (in the order of the blocks of rows and columns of A_TF). The
		//   bases must be hierarchical (ordered by degree: for each d, the first functions span P^d, or Q^d whose degree is
		//   the largest partial degree; not needed with PCoarsening::ToZero) with a constant first function phi_0:
		//   1 = c phi_0, i.e. c = 1/phi_0. phi_0 may differ from one cell or face to another (e.g. L2-orthonormal bases:
		//   c = sqrt(|T|) on a cell T); c = 1 everywhere if phi_0 = 1.

		//----------------------------- Setups and solves -----------------------------//

		// Setup for a symmetric positive definite system, from its blocks (A_FT = A_TF^T), without the DoFs of the
		// Dirichlet faces. Computes the condensed matrix A = A_FF - A_FT A_TT^{-1} A_TF (CondensedMatrix()), and keeps what
		// Solve(b_T, b_F, x_T, x_F) needs: the blocks may be freed after the setup.
		// A_TT (cells x cells), A_TF (cells x faces), A_FF (faces x faces): the DoFs of each cell (resp. face) are
		// contiguous, in the order of the rows (resp. columns) of A_TF. A_TT and A_FF are read from their lower triangular
		// part (full or lower-only matrices). The blocks of A_TT (one per cell) must be positive definite: they are
		// factorized by Cholesky (fhhos4::Error otherwise).
		// On several processes (Comm), each one gives the rows of its own cells and faces: A_TT and A_TF the rows of its
		// cells, A_FF those of its faces, in full (both triangles); the columns of A_TF and A_FF are numbered locally,
		// in one numbering of the faces of the process: its own faces (those of the rows of A_FF, in their order), then
		// its ghosts, the other faces of its rows (of its cells, and those coupled to its faces in A_FF through the cells
		// of other processes), whose global numbers are faceGhosts; the interpolations of 1 on its cells and faces. Each
		// face belongs to a process that owns one of its cells (e.g. the one of the smaller rank); a process may own no
		// cell nor face. Then Solve() with the vectors of the cells and faces of the process. On one process, faceGhosts
		// is empty.
		template <typename TT, typename TF, typename FF>
		void SetupSymmetric(const Eigen::SparseCompressedBase<TT>& A_TT, const Eigen::SparseCompressedBase<TF>& A_TF, const Eigen::SparseCompressedBase<FF>& A_FF,
		                    const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne,
		                    const GhostNumbers& faceGhosts = GhostNumbers())
		{
			Setup(true, View(A_TT), View(A_TF), nullptr, View(A_FF), In(cellInterpOfOne), In(faceInterpOfOne), faceGhosts);
		}

		// Setup for a non-symmetric system (e.g. convection-diffusion), whose symmetric part is positive definite: as
		// SetupSymmetric(), with the block A_FT (faces x cells) and all the blocks read in full: the same
		// algorithm runs on the matrix itself. The couplings of the aggregation are computed from the averaged traces of
		// the blocks of A_TF and A_FT; the coarse operators are the Galerkin products P^T A P (as AGMG: Y. Notay, "An
		// aggregation-based algebraic multigrid method", ETNA 37, 2010); the Gauss-Seidel smoothers follow the downwind
		// order of each level (convection), in the pre- and post-smoothing; the K-cycle iterates with GCR. Krylov must not
		// be FCG, nor CoarseSolver Cholesky, nor the smoothers Chebyshev.
		// The solver keeps a reference to A_FT, read by each Solve(b_T, b_F, x_T, x_F) (b_F - A_FT A_TT^{-1} b_T): A_FT must
		// stay alive and unchanged until the next setup or the destruction of the solver (not a temporary matrix). The
		// other blocks may be freed after the setup.
		// On several processes (Comm), as SetupSymmetric(), A_TT and A_FF read in full, and A_FT the rows of the faces of
		// the process, its columns numbered locally: the cells of the process (those of the rows of A_TT, in their order),
		// then its ghost cells, the other cells of its rows (of other processes), whose global numbers are cellGhosts (as
		// faceGhosts: the cells of the process of rank 0 first, etc.). On one process, faceGhosts and cellGhosts are empty.
		// See README_library.md.
		template <typename TT, typename TF, typename FT, typename FF>
		void SetupNonsymmetric(const Eigen::SparseCompressedBase<TT>& A_TT, const Eigen::SparseCompressedBase<TF>& A_TF, const Eigen::SparseCompressedBase<FT>& A_FT,
		                       const Eigen::SparseCompressedBase<FF>& A_FF,
		                       const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne,
		                       const GhostNumbers& faceGhosts = GhostNumbers(), const GhostNumbers& cellGhosts = GhostNumbers())
		{
			detail::CSRView a_FT = View(A_FT);
			Setup(false, View(A_TT), View(A_TF), &a_FT, View(A_FF), In(cellInterpOfOne), In(faceInterpOfOne), faceGhosts, cellGhosts);
		}

		template <typename TT, typename TF, typename FF, int Options, typename Index>
		void SetupNonsymmetric(const Eigen::SparseCompressedBase<TT>&, const Eigen::SparseCompressedBase<TF>&, Eigen::SparseMatrix<double, Options, Index>&& A_FT,
		                       const Eigen::SparseCompressedBase<FF>&, const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&,
		                       const GhostNumbers& = GhostNumbers(), const GhostNumbers& = GhostNumbers()) = delete;

		// Solves the system, after SetupSymmetric() or SetupNonsymmetric(): condenses the right-hand
		// side, b = b_F - A_FT A_TT^{-1} b_T, solves A x_F = b, then computes x_T = A_TT^{-1} (b_T - A_TF x_F).
		// x_F: the initial guess on input (e.g. zero), the solution on output; x_T: output. The result is that of
		// A x_F = b. Does not throw if the solver has not converged after MaxIterations iterations: see the result.
		// On several processes: the vectors of the cells and faces of the process; the result, of the whole system.
		Result Solve(const Eigen::Ref<const Eigen::VectorXd>& b_T, const Eigen::Ref<const Eigen::VectorXd>& b_F,
		             Eigen::Ref<Eigen::VectorXd> x_T, Eigen::Ref<Eigen::VectorXd> x_F)
		{
			return SolveUncondensed(In(b_T), In(b_F), Out(x_T), Out(x_F));
		}

		// Solves the condensed system A x = b (A = A_FF - A_FT A_TT^{-1} A_TF: CondensedMatrix() on one process), after a
		// setup, e.g. for right-hand sides condensed by the calling code. x: the initial guess on input (e.g. zero), the solution on
		// output. Does not throw if the solver has not converged after MaxIterations iterations: see the result.
		Result Solve(const Eigen::Ref<const Eigen::VectorXd>& b, Eigen::Ref<Eigen::VectorXd> x)
		{
			return SolveCondensed(In(b), Out(x));
		}

		// Solves A x = b from x = 0
		Eigen::VectorXd Solve(const Eigen::Ref<const Eigen::VectorXd>& b)
		{
			Eigen::VectorXd x = Eigen::VectorXd::Zero(b.rows());
			Solve(b, x);
			return x;
		}

		//------------------------------------------------------------------//

		// The condensed matrix A of the last setup, both triangles, computed by SetupSymmetric() or SetupNonsymmetric().
		// Valid until the next setup or the destruction of the solver. On one process only.
		SparseMatrixMap CondensedMatrix() const
		{
			detail::CSRView A = CondensedMatrixView();
			return SparseMatrixMap(A.Rows, A.Cols, A.NonZeros(), A.RowPointers, A.ColumnIndices, A.Values);
		}

		bool IsSetUp() const;
		int NumberOfLevels() const;

	protected:
		// The implementation of the setups, on the arrays of the matrices (also for the Python module). A_FT: null for a
		// symmetric matrix; faceGhosts, cellGhosts (A_FT): on several processes.
		void Setup(bool symmetric, const detail::CSRView& A_TT, const detail::CSRView& A_TF, const detail::CSRView* A_FT, const detail::CSRView& A_FF,
		           detail::ConstVectorView cellInterpOfOne, detail::ConstVectorView faceInterpOfOne, const GhostNumbers& faceGhosts,
		           const GhostNumbers& cellGhosts = GhostNumbers());
		detail::CSRView CondensedMatrixView() const;

		static detail::ConstVectorView In(const Eigen::Ref<const Eigen::VectorXd>& v)
		{
			return { v.data(), (std::int64_t)v.size() };
		}

		static detail::VectorView Out(Eigen::Ref<Eigen::VectorXd> v)
		{
			return { v.data(), (std::int64_t)v.size() };
		}

	private:
		struct Impl;
		std::unique_ptr<Impl> _impl;

		// The arrays of a matrix given to a setup: an Eigen sparse matrix in the format of the library, read in place
		template <typename Derived>
		static detail::CSRView View(const Eigen::SparseCompressedBase<Derived>& M)
		{
			constexpr bool format = Derived::IsRowMajor && std::is_same<typename Derived::StorageIndex, SparseIndex>::value;
			static_assert(format, "fhhos4: the matrices are given in the format of the library: fhhos4::SparseMatrix (row-major, with the index type fhhos4::SparseIndex), or an Eigen::Map of one. Convert this one by Eigen (a copy): fhhos4::SparseMatrix B = A; or, if it is symmetric and column-major, view it with fhhos4::TransposedView(A) (no copy).");
			static_assert(std::is_same<typename Derived::Scalar, double>::value, "fhhos4: matrices of doubles only");
			if constexpr (format)
				return detail::View<SparseIndex>(M);
			else
				return detail::CSRView();
		}

		void SetupDistributed(bool symmetric, const detail::CSRView& A_TT, const detail::CSRView& A_TF, const detail::CSRView* A_FT, const detail::CSRView& A_FF,
		                      detail::ConstVectorView cellInterpOfOne, detail::ConstVectorView faceInterpOfOne, const GhostNumbers& faceGhosts,
		                      const GhostNumbers& cellGhosts);
		Result SolveCondensed(detail::ConstVectorView b, detail::VectorView x);
		Result SolveUncondensed(detail::ConstVectorView b_T, detail::ConstVectorView b_F, detail::VectorView x_T, detail::VectorView x_F);
	};

	//------------------------------------------------------------------------------------------------------------------//
	//                       Helpers: the inputs of the setups from other formats (copies)                            //
	//------------------------------------------------------------------------------------------------------------------//

	// The DoFs of each cell (or face) in the numbering of a global matrix, read in place: those of the i-th are
	// DoFs[Offsets[i]], ..., DoFs[Offsets[i + 1] - 1], in the order of its basis. Index: std::int32_t or std::int64_t.
	template <typename Index>
	struct DoFLists
	{
		static_assert(std::is_same<Index, std::int32_t>::value || std::is_same<Index, std::int64_t>::value, "fhhos4::DoFLists: the indices are std::int32_t or std::int64_t");

		std::int64_t Count = 0;         // number of cells (or faces)
		const Index* Offsets = nullptr; // Count + 1, Offsets[0] = 0
		const Index* DoFs = nullptr;    // Offsets[Count]
	};

	// The blocks of a hybrid system: see ExtractSymmetricBlocks()
	struct HybridBlocks
	{
		SparseMatrix A_TT, A_TF, A_FT, A_FF; // A_FT: empty (0 x 0) for a symmetric system
	};

	namespace detail
	{
		// The blocks in the arrays of the library
		struct HybridBlocks
		{
			CSRMatrix A_TT, A_TF, A_FT, A_FF;
		};

		template <typename RowPointer, typename Column, typename Index>
		FHHOS4_API HybridBlocks ExtractSymmetricBlocks(const BasicCSRView<RowPointer, Column>& M, const DoFLists<Index>& cells, const DoFLists<Index>& faces);
		template <typename RowPointer, typename Column, typename Index>
		FHHOS4_API HybridBlocks ExtractNonsymmetricBlocks(const BasicCSRView<RowPointer, Column>& M, const DoFLists<Index>& cells, const DoFLists<Index>& faces);
		template <typename Index>
		FHHOS4_API void Gather(ConstVectorView x, const DoFLists<Index>& cells, const DoFLists<Index>& faces, VectorView x_T, VectorView x_F);
		template <typename Index>
		FHHOS4_API void Scatter(ConstVectorView x_T, ConstVectorView x_F, const DoFLists<Index>& cells, const DoFLists<Index>& faces, VectorView x);
		// A copy: M in the CSR format of the library, from other index types
		template <typename RowPointer, typename Column>
		FHHOS4_API CSRMatrix ToCSR(const BasicCSRView<RowPointer, Column>& M);
		// A copy: M^T in the CSR format of the library
		template <typename Index>
		FHHOS4_API CSRMatrix Transpose(const BasicCSRView<Index, Index>& M);

		// The blocks copied into Eigen matrices
		inline fhhos4::HybridBlocks ToEigen(const HybridBlocks& blocks)
		{
			return { ToEigen(blocks.A_TT), ToEigen(blocks.A_TF), ToEigen(blocks.A_FT), ToEigen(blocks.A_FF) };
		}

		// The arrays of the global matrix M of the helpers: compressed, row-major, of doubles, with the indices
		// std::int32_t or std::int64_t
		template <typename Derived>
		BasicCSRView<typename Derived::StorageIndex, typename Derived::StorageIndex> GlobalMatrixView(const Eigen::SparseCompressedBase<Derived>& M)
		{
			return View<typename Derived::StorageIndex>(M);
		}
	}

	// A copy: the blocks of the uncondensed system, extracted from its global matrix M (square), for SetupSymmetric():
	// A_TT, A_TF, A_FF, in one parallel pass (two on the rows of M). M: an Eigen sparse matrix, compressed, row-major, with
	// the index type std::int32_t or std::int64_t (e.g. an Eigen::Map of the arrays of another code). The cells (resp.
	// faces) in the order of their lists, the DoFs of each in the order of its list (e.g. reordered by degree); the DoFs
	// in no list are dropped (e.g. those of the Dirichlet faces: move their contribution to the right-hand side first).
	// Each DoF in one list at most. Symmetric matrix: A_TT and A_FF lower triangular (what SetupSymmetric() reads).
	template <typename Derived, typename Index>
	HybridBlocks ExtractSymmetricBlocks(const Eigen::SparseCompressedBase<Derived>& M, const DoFLists<Index>& cells, const DoFLists<Index>& faces)
	{
		static_assert(Derived::IsRowMajor, "fhhos4::ExtractSymmetricBlocks(): a row-major matrix (convert a column-major one by Eigen, a copy)");
		static_assert(std::is_same<typename Derived::Scalar, double>::value, "fhhos4: matrices of doubles only");
		return detail::ToEigen(detail::ExtractSymmetricBlocks(detail::GlobalMatrixView(M), cells, faces));
	}

	// A copy: as ExtractSymmetricBlocks(), the four blocks in full, for SetupNonsymmetric()
	template <typename Derived, typename Index>
	HybridBlocks ExtractNonsymmetricBlocks(const Eigen::SparseCompressedBase<Derived>& M, const DoFLists<Index>& cells, const DoFLists<Index>& faces)
	{
		static_assert(Derived::IsRowMajor, "fhhos4::ExtractNonsymmetricBlocks(): a row-major matrix (convert a column-major one by Eigen, a copy)");
		static_assert(std::is_same<typename Derived::Scalar, double>::value, "fhhos4: matrices of doubles only");
		return detail::ToEigen(detail::ExtractNonsymmetricBlocks(detail::GlobalMatrixView(M), cells, faces));
	}

	// A copy: the vectors of the cells and faces DoFs (e.g. b_T and b_F) from a vector x in the numbering of the global
	// matrix (b), with the lists of ExtractSymmetricBlocks()
	template <typename Index>
	void Gather(const Eigen::Ref<const Eigen::VectorXd>& x, const DoFLists<Index>& cells, const DoFLists<Index>& faces,
	            Eigen::Ref<Eigen::VectorXd> x_T, Eigen::Ref<Eigen::VectorXd> x_F)
	{
		detail::Gather<Index>({ x.data(), (std::int64_t)x.size() }, cells, faces, { x_T.data(), (std::int64_t)x_T.size() }, { x_F.data(), (std::int64_t)x_F.size() });
	}

	// A copy: the other way around, x_T and x_F into the vector x in the numbering of the global matrix (its DoFs in no
	// list unchanged, e.g. those of the Dirichlet faces)
	template <typename Index>
	void Scatter(const Eigen::Ref<const Eigen::VectorXd>& x_T, const Eigen::Ref<const Eigen::VectorXd>& x_F, const DoFLists<Index>& cells, const DoFLists<Index>& faces,
	             Eigen::Ref<Eigen::VectorXd> x)
	{
		detail::Scatter<Index>({ x_T.data(), (std::int64_t)x_T.size() }, { x_F.data(), (std::int64_t)x_F.size() }, cells, faces, { x.data(), (std::int64_t)x.size() });
	}

	// No copy: A^T in the format of the library, a map of the arrays of a column-major matrix A (compressed, with the
	// index type SparseIndex), which is A itself if A is symmetric, up to the rounding errors of its assembly (its
	// asymmetry): for the symmetric inputs of the setups (A_TT, A_FF, and A of a symmetric system) stored column-major.
	// Not of a temporary.
	template <typename Derived>
	SparseMatrixMap TransposedView(const Eigen::SparseCompressedBase<Derived>& A)
	{
		constexpr bool format = !Derived::IsRowMajor && std::is_same<typename Derived::StorageIndex, SparseIndex>::value;
		static_assert(!Derived::IsRowMajor, "fhhos4::TransposedView(): for a column-major matrix (give a row-major one in the format of the library as it is)");
		static_assert(std::is_same<typename Derived::StorageIndex, SparseIndex>::value, "fhhos4::TransposedView(): the index type of the matrix is not fhhos4::SparseIndex: convert it by Eigen (a copy): fhhos4::SparseMatrix B = A;");
		static_assert(std::is_same<typename Derived::Scalar, double>::value, "fhhos4: matrices of doubles only");
		if (!A.isCompressed())
			throw Error("fhhos4::TransposedView(): the Eigen matrix is not compressed: call its makeCompressed() first.");
		if constexpr (format)
			return SparseMatrixMap(A.cols(), A.rows(), A.nonZeros(), A.outerIndexPtr(), A.innerIndexPtr(), A.valuePtr());
		else
			return SparseMatrixMap(0, 0, 0, nullptr, nullptr, nullptr);
	}

	template <typename Derived>
	SparseMatrixMap TransposedView(const Eigen::SparseCompressedBase<Derived>&&) = delete;
}
