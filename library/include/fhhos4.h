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
// Solves a hybrid system (cell and face unknowns), symmetric positive definite (SetupSymmetric...()) or not
// (SetupNonsymmetric...(), e.g. convection-diffusion),
//     [ A_TT  A_TF ] [x_T]   [b_T]
//     [ A_FT  A_FF ] [x_F] = [b_F],
// by static condensation: the AMG solves the condensed system A x_F = b, A = A_FF - A_FT A_TT^{-1} A_TF and
// b = b_F - A_FT A_TT^{-1} b_T. Two interfaces:
// - the uncondensed system, the interface of choice: SetupSymmetric() or SetupNonsymmetric() with the blocks, then
//   Solve(b_T, b_F, x_T, x_F). The library condenses the system and computes the cell unknowns: no static condensation
//   in the calling code;
// - the condensed system, to integrate the library with the least effort into a code that condenses statically
//   already: Setup...Condensed() with A and the blocks, then Solve(b, x).
// See README_library.md.
//
// The Krylov method preconditioned by the AMG is part of the library (Solver::Krylov), which does not give the multigrid
// cycle alone: with the default K-cycle, the AMG is a non-linear preconditioner (it solves the coarse levels with inner
// FCG or GCR iterations), which requires a flexible Krylov method: the flexible conjugate gradient (FCG, the default for
// a symmetric matrix) or GCR.
//
// The inputs are read in place, never copied nor converted: the sparse matrices in the CSR format of the library
// (CSRView: row-major, with the index type SparseIndex), the vectors as arrays of doubles. The helpers at the end of
// this file convert the other formats (other index types, column-major matrices, the global matrix of the system
// instead of its blocks): they copy them. Assembling the matrices in the format of the library avoids the copies.
//
// The functions compiled in the shared library fhhos4 take arrays only (views): the Eigen objects of this header (the
// Eigen matrices converted to views, the Eigen::Ref of the vectors) are those of the calling code, compiled with its
// own version of Eigen. The library never allocates memory that the calling code frees, nor the other way around.
//
// Compiled with MPI (FHHOS4_MPI, -DENABLE_MPI=ON), a Solver works on the processes of its communicator (Solver::Comm),
// on a system distributed by rows: each process gives the rows of its own cells and faces, the columns of the faces
// numbered locally, its own faces first, then the other faces of its cells (its ghosts), whose global numbers it gives
// too (GhostNumbers). See the setups below, and README_library.md.
namespace fhhos4
{
	// The index type of the sparse matrices of the library: int (32 bits), or 64 bits if fhhos4 is compiled with
	// -DSMALL_INDEX=OFF
	using SparseIndex = FHHOS4_SPARSE_INDEX;

	// The sparse matrices of the library as Eigen matrices: row-major, with the index type SparseIndex. Compressed, they
	// convert to a CSRView without copy.
	using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, SparseIndex>;

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

	// A sparse matrix in the CSR format (compressed sparse row), not owned: the arrays of the calling code, read in place.
	// The non-zeros of the row i are (ColumnIndices[p], Values[p]) for RowPointers[i] <= p < RowPointers[i + 1], with
	// RowPointers[0] = 0 and the columns strictly increasing in each row (no duplicates). RowPointer and Column:
	// std::int32_t or std::int64_t, independently; the setups take a CSRView, the index type of the library for both.
	template <typename RowPointer, typename Column>
	struct BasicCSRView
	{
		static_assert((std::is_same<RowPointer, std::int32_t>::value || std::is_same<RowPointer, std::int64_t>::value) &&
		              (std::is_same<Column, std::int32_t>::value || std::is_same<Column, std::int64_t>::value),
		              "fhhos4::BasicCSRView: the indices are std::int32_t or std::int64_t");

		std::int64_t Rows = 0;
		std::int64_t Cols = 0;
		const RowPointer* RowPointers = nullptr; // Rows + 1
		const Column* ColumnIndices = nullptr;   // RowPointers[Rows]
		const double* Values = nullptr;          // RowPointers[Rows]

		BasicCSRView() = default;

		BasicCSRView(std::int64_t rows, std::int64_t cols, const RowPointer* rowPointers, const Column* columnIndices, const double* values)
			: Rows(rows), Cols(cols), RowPointers(rowPointers), ColumnIndices(columnIndices), Values(values)
		{}

		// The arrays of an Eigen sparse matrix (Eigen::SparseMatrix, or an Eigen::Map of one): compressed, row-major, with
		// these index types, e.g. fhhos4::SparseMatrix for a CSRView. A column-major matrix stores its transpose: see
		// ToCSR() (a copy) and TransposedView() (no copy, for a symmetric matrix).
		template <typename Derived>
		BasicCSRView(const Eigen::SparseCompressedBase<Derived>& M)
		{
			constexpr bool rowMajor = Derived::IsRowMajor;
			constexpr bool sameIndices = std::is_same<typename Derived::StorageIndex, RowPointer>::value && std::is_same<typename Derived::StorageIndex, Column>::value;
			static_assert(rowMajor, "fhhos4: a column-major matrix stores its transpose. Convert it with fhhos4::ToCSR() (a copy), or, if it is symmetric, view it with fhhos4::TransposedView() (no copy).");
			static_assert(sameIndices, "fhhos4: the index type of the matrix is not that of the view (fhhos4::SparseIndex for a CSRView). Use fhhos4::SparseMatrix, or convert the matrix with fhhos4::ToCSR() (a copy).");
			static_assert(std::is_same<typename Derived::Scalar, double>::value, "fhhos4: matrices of doubles only");
			if constexpr (rowMajor && sameIndices)
			{
				if (!M.isCompressed())
					throw Error("fhhos4: the Eigen matrix is not compressed: call its makeCompressed() first.");
				Rows = M.rows();
				Cols = M.cols();
				RowPointers = M.outerIndexPtr();
				ColumnIndices = M.innerIndexPtr();
				Values = M.valuePtr();
			}
		}

		// Not of a temporary matrix: the view would outlive it
		template <typename Derived>
		BasicCSRView(const Eigen::SparseCompressedBase<Derived>&&) = delete;

		std::int64_t NonZeros() const
		{
			return RowPointers ? (std::int64_t)RowPointers[Rows] : 0;
		}

		// The view as an Eigen matrix (RowPointer and Column of the same type), e.g. CondensedMatrix().EigenMap()
		Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor, RowPointer>> EigenMap() const
		{
			static_assert(std::is_same<RowPointer, Column>::value, "fhhos4: EigenMap() needs row pointers and columns of the same type");
			return Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor, RowPointer>>(Rows, Cols, NonZeros(), RowPointers, ColumnIndices, Values);
		}
	};

	// The matrices of the setups: the index type of the library
	using CSRView = BasicCSRView<SparseIndex, SparseIndex>;

	// A sparse matrix in the CSR format of the library, owned by the library (allocated and freed by it): the result of
	// the helpers that copy (ToCSR(), ExtractSymmetricBlocks()...). Converts to a CSRView, valid while the object lives
	// (not of a temporary).
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
		operator CSRView() const& { return View(); }
		operator CSRView() const&& = delete;

	private:
		std::unique_ptr<Impl> _impl;
	};

	// On several processes (MPI): the global numbers of the ghosts of a process (the faces of its rows owned by other
	// processes), in the order of their local numbers (after those of its own faces), read in place. The global numbers:
	// the faces of the process of rank 0 first, then those of the process of rank 1, etc., each process's in the order
	// of its local numbers (so that the faces of each process are numbered contiguously). Any order of the ghosts, no
	// repetition.
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
		BiCGSTAB, // Eigen's BiCGSTAB. Not flexible: it needs a linear preconditioner, a V- or W-cycle (not K). Not on
		          // several processes (MPI)
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

	// The p-coarsening (Solver::PCoarsening): the degrees of the p-levels, from FaceDegree down to 0, before the h-levels
	// (aggregation) of degree 0. A p-level of degree d keeps the first functions of each block of the finer level: those
	// of the polynomials of degree <= d if the bases are ordered by degree.
	enum class PCoarsening
	{
		Minus2,    // k, k - 2, ..., 0 (fhhos4 -p-cs -2)
		Minus1,    // k, k - 1, ..., 0 (-p-cs -1)
		DivideBy2, // k, k / 2, ..., 0 (-p-cs /2)
		ToZero     // k, 0 (-p-cs =0): the only p-level keeps the first function of each block, the constant one: for bases
		           // not ordered by degree. Usually more iterations than the others at k >= 3
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

	// The solver of the coarsest level (Solver::CoarseSolver)
	enum class CoarseSolver
	{
		Default,  // Cholesky for a symmetric matrix (SetupSymmetric...()), LU otherwise (SetupNonsymmetric...())
		Cholesky, // sparse Cholesky factorization, for a symmetric positive definite matrix
		LU        // sparse LU factorization
	};

	struct Result
	{
		int Iterations = 0;
		double RelativeResidual = 0; // |b - A x| / |b|, in the condensed system
		bool Converged = false;      // RelativeResidual < Tolerance
	};

	// Set the parameters, then call a setup function, then Solve(). The parameters are read by the setup (Tolerance and
	// MaxIterations by Solve()). The defaults are those of the program fhhos4 (-s fcguamg), used in the papers, except the
	// smoothers (PreSmoother, and the number of its iterations).
	// Not thread-safe: one object must not be used by several threads at the same time.
	class FHHOS4_API Solver
	{
	public:
		//------------- Discretization: to set -------------//

		int Dimension = -1;  // 2 or 3
		int FaceDegree = -1; // k: degree of the polynomials on the faces
		int CellDegree = -1; // l: degree of the polynomials in the cells (k - 1, k or k + 1 in HHO)

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
		double CoarseningFactor = 3.8;     // target ratio of the numbers of cells of two successive h-levels
		int CoarseMatrixMaxSize = 1000;    // the coarsening stops when the matrix has at most this many rows
		// Solver of the coarsest level. Default: Cholesky for a symmetric matrix, LU otherwise.
		fhhos4::CoarseSolver CoarseSolver = fhhos4::CoarseSolver::Default;

		//-------------------- Execution -------------------//

		int Threads = 0;   // number of OpenMP threads during the calls of the library; 0: the current OpenMP setting, but
		                   // at most one thread per physical core (no hyper-threading) unless OMP_NUM_THREADS or OpenMP's
		                   // binding (OMP_PROC_BIND, OMP_PLACES) is set. On several processes (MPI): the cores divided
		                   // among the processes of Comm that run on the same CPUs (not bound to distinct cores)
		int Verbosity = 0; // 0: nothing; 1: one line per setup (number of levels, time) and per solve (iterations,
		                   // residual, time); 2: all the messages of fhhos4, with the residual at each iteration.
		                   // On several processes: the messages of the process of rank 0.

#if FHHOS4_MPI
		// The processes of the solver (the library is compiled with MPI): to set, no default. MPI_COMM_SELF: the system of
		// each process, solved by it alone, as without MPI. A communicator of several processes: one system distributed by
		// rows over them (see SetupSymmetric()); every setup and solve is then collective, called by all the processes of
		// the communicator, with the same parameters.
		MPI_Comm Comm = MPI_COMM_NULL;
#endif

		Solver();
		~Solver();
		Solver(Solver&&) noexcept;
		Solver& operator=(Solver&&) noexcept;
		Solver(const Solver&) = delete;
		Solver& operator=(const Solver&) = delete;

		// The inputs of the setups:
		// - the matrices (CSRView): read in place, in the CSR format of the library, never copied. An Eigen matrix
		//   (compressed, row-major, with the index type SparseIndex: fhhos4::SparseMatrix) converts to a CSRView, a
		//   CSRMatrix too; other formats: see the helpers at the end of this file (they copy). Not temporaries. The setups
		//   check the sizes and structure of the inputs against Dimension and the degrees, and the CSR format (one pass on
		//   the indices);
		// - cellInterpOfOne, faceInterpOfOne: the interpolation of the function 1 on the polynomial bases of the cells and
		//   faces, one coefficient per cell and per face (in the order of the blocks of rows and columns of A_TF). The
		//   bases must be hierarchical (ordered by degree; not needed with PCoarsening::ToZero) with a constant first
		//   function phi_0: 1 = c phi_0, i.e. c = 1/phi_0. phi_0 may differ from one cell or face to another (e.g.
		//   L2-orthonormal bases: c = sqrt(|T|) on a cell T); c = 1 everywhere if phi_0 = 1.

		//----------- Uncondensed system: the interface of choice -----------//

		// Setup for a symmetric positive definite system, from its blocks (A_FT = A_TF^T), without the DoFs of the
		// Dirichlet faces. Computes the condensed matrix A = A_FF - A_FT A_TT^{-1} A_TF (CondensedMatrix()), and keeps what
		// Solve(b_T, b_F, x_T, x_F) needs: the blocks may be freed after the setup.
		// A_TT (cells x cells), A_TF (cells x faces), A_FF (faces x faces): the DoFs of each cell (resp. face) are
		// contiguous, in the order of the rows (resp. columns) of A_TF. A_TT and A_FF are read from their lower triangular
		// part (full or lower-only matrices).
		// On several processes (Comm), each one gives the rows of its own cells and faces: A_TT and A_TF the rows of its
		// cells, A_FF those of its faces, in full (both triangles); the columns of A_TF and A_FF are numbered locally,
		// in one numbering of the faces of the process: its own faces (those of the rows of A_FF, in their order), then
		// its ghosts, the other faces of its rows (of its cells, and those coupled to its faces in A_FF through the cells
		// of other processes), whose global numbers are faceGhosts; the interpolations of 1 on its cells and faces. Each
		// face belongs to a process that owns one of its cells (e.g. the one of the smaller rank); a process may own no
		// cell nor face. Then Solve() with the vectors of the cells and faces of the process. On one process, faceGhosts
		// is empty.
		void SetupSymmetric(const CSRView& A_TT, const CSRView& A_TF, const CSRView& A_FF,
		                    const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne,
		                    const GhostNumbers& faceGhosts = GhostNumbers())
		{
			Setup(true, nullptr, A_TT, A_TF, nullptr, &A_FF, In(cellInterpOfOne), In(faceInterpOfOne), faceGhosts);
		}

		// Setup for a non-symmetric system (e.g. convection-diffusion), whose symmetric part is positive definite: as
		// SetupSymmetric() (on one process only, for now), with the block A_FT (faces x cells) and all the blocks read in
		// full: the same
		// algorithm runs on the matrix itself. The couplings of the aggregation are computed from the averaged traces of
		// the blocks of A_TF and A_FT; the coarse operators are the Galerkin products P^T A P (as AGMG: Y. Notay, "An
		// aggregation-based algebraic multigrid method", ETNA 37, 2010); the Gauss-Seidel smoothers follow the downwind
		// order of each level (convection), in the pre- and post-smoothing; the K-cycle iterates with GCR. Krylov must not
		// be FCG, nor CoarseSolver Cholesky, nor the smoothers Chebyshev.
		// The solver keeps a reference to A_FT, read by each Solve(b_T, b_F, x_T, x_F) (b_F - A_FT A_TT^{-1} b_T): A_FT must
		// stay alive and unchanged until the next setup or the destruction of the solver. The other blocks may be freed
		// after the setup.
		// See README_library.md.
		void SetupNonsymmetric(const CSRView& A_TT, const CSRView& A_TF, const CSRView& A_FT, const CSRView& A_FF,
		                       const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
		{
			Setup(false, nullptr, A_TT, A_TF, &A_FT, &A_FF, In(cellInterpOfOne), In(faceInterpOfOne), GhostNumbers());
		}

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

		//---- Condensed system: to integrate the library into a code that condenses statically ----//
		// On one process only (on several, give the uncondensed system: SetupSymmetric()).

		// Setup for a symmetric positive definite matrix, from the condensed matrix and the blocks: as
		// SetupSymmetric(), with A (faces x faces), the condensed matrix A_FF - A_FT A_TT^{-1} A_TF, both
		// triangles, and A_FF optional (the algorithm only uses A_TT and A_TF, whose rows each belong to one cell).
		// The solver keeps a reference to A: A must stay alive and unchanged until the next setup or the destruction of
		// the solver. The blocks are only read by the setup: they may be freed after it.
		void SetupSymmetricCondensed(const CSRView& A, const CSRView& A_TT, const CSRView& A_TF, const CSRView& A_FF,
		                             const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
		{
			Setup(true, &A, A_TT, A_TF, nullptr, &A_FF, In(cellInterpOfOne), In(faceInterpOfOne), GhostNumbers());
		}

		void SetupSymmetricCondensed(const CSRView& A, const CSRView& A_TT, const CSRView& A_TF,
		                             const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
		{
			Setup(true, &A, A_TT, A_TF, nullptr, nullptr, In(cellInterpOfOne), In(faceInterpOfOne), GhostNumbers());
		}

		// Setup for a non-symmetric matrix: as SetupNonsymmetric(), with the condensed matrix A (see
		// SetupSymmetricCondensed()) and A_FF optional. A_FT is only read by the setup, as the other blocks.
		void SetupNonsymmetricCondensed(const CSRView& A, const CSRView& A_TT, const CSRView& A_TF, const CSRView& A_FT, const CSRView& A_FF,
		                                const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
		{
			Setup(false, &A, A_TT, A_TF, &A_FT, &A_FF, In(cellInterpOfOne), In(faceInterpOfOne), GhostNumbers());
		}

		void SetupNonsymmetricCondensed(const CSRView& A, const CSRView& A_TT, const CSRView& A_TF, const CSRView& A_FT,
		                                const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
		{
			Setup(false, &A, A_TT, A_TF, &A_FT, nullptr, In(cellInterpOfOne), In(faceInterpOfOne), GhostNumbers());
		}

		// Solves the condensed system A x = b, after any setup. x: the initial guess on input (e.g. zero), the solution on
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

		// The condensed matrix A of the last setup, both triangles: computed by SetupSymmetric() and SetupNonsymmetric(),
		// or the one given to SetupSymmetricCondensed() and SetupNonsymmetricCondensed(). Valid until the next setup or the
		// destruction of the solver. As an Eigen matrix: CondensedMatrix().EigenMap(). On one process only.
		CSRView CondensedMatrix() const;

		bool IsSetUp() const;
		int NumberOfLevels() const;

	private:
		struct Impl;
		std::unique_ptr<Impl> _impl;

		static detail::ConstVectorView In(const Eigen::Ref<const Eigen::VectorXd>& v)
		{
			return { v.data(), (std::int64_t)v.size() };
		}

		static detail::VectorView Out(Eigen::Ref<Eigen::VectorXd> v)
		{
			return { v.data(), (std::int64_t)v.size() };
		}

		// The setups. A: null for SetupSymmetric() and SetupNonsymmetric(); A_FT: null for a symmetric matrix; A_FF: null
		// if not given; faceGhosts: on several processes.
		void Setup(bool symmetric, const CSRView* A, const CSRView& A_TT, const CSRView& A_TF, const CSRView* A_FT, const CSRView* A_FF,
		           detail::ConstVectorView cellInterpOfOne, detail::ConstVectorView faceInterpOfOne, const GhostNumbers& faceGhosts);
		void SetupDistributed(bool symmetric, const CSRView* A, const CSRView& A_TT, const CSRView& A_TF, const CSRView* A_FT, const CSRView* A_FF,
		                      detail::ConstVectorView cellInterpOfOne, detail::ConstVectorView faceInterpOfOne, const GhostNumbers& faceGhosts);
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

	// The blocks of a hybrid system, owned by the library: see ExtractSymmetricBlocks()
	struct HybridBlocks
	{
		CSRMatrix A_TT, A_TF, A_FT, A_FF; // A_FT: empty (0 x 0) for a symmetric system
	};

	// A copy: the blocks of the uncondensed system, extracted from its global matrix M (square, in the CSR format, of any
	// index types), for SetupSymmetric(): A_TT, A_TF, A_FF, in one parallel pass (two on the rows of M). The cells (resp.
	// faces) in the order of their lists, the DoFs of each in the order of its list (e.g. reordered by degree); the DoFs
	// in no list are dropped (e.g. those of the Dirichlet faces: move their contribution to the right-hand side first).
	// Each DoF in one list at most. Symmetric matrix: A_TT and A_FF lower triangular (what SetupSymmetric() reads).
	template <typename RowPointer, typename Column, typename Index>
	FHHOS4_API HybridBlocks ExtractSymmetricBlocks(const BasicCSRView<RowPointer, Column>& M, const DoFLists<Index>& cells, const DoFLists<Index>& faces);

	// A copy: as ExtractSymmetricBlocks(), the four blocks in full, for SetupNonsymmetric()
	template <typename RowPointer, typename Column, typename Index>
	FHHOS4_API HybridBlocks ExtractNonsymmetricBlocks(const BasicCSRView<RowPointer, Column>& M, const DoFLists<Index>& cells, const DoFLists<Index>& faces);

	namespace detail
	{
		template <typename Index>
		FHHOS4_API void Gather(ConstVectorView x, const DoFLists<Index>& cells, const DoFLists<Index>& faces, VectorView x_T, VectorView x_F);
		template <typename Index>
		FHHOS4_API void Scatter(ConstVectorView x_T, ConstVectorView x_F, const DoFLists<Index>& cells, const DoFLists<Index>& faces, VectorView x);
		// A copy: M^T in the CSR format of the library
		template <typename Index>
		FHHOS4_API CSRMatrix Transpose(const BasicCSRView<Index, Index>& M);
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

	// A copy: M in the CSR format of the library, from other index types
	template <typename RowPointer, typename Column>
	FHHOS4_API CSRMatrix ToCSR(const BasicCSRView<RowPointer, Column>& M);

	// A copy: an Eigen sparse matrix (compressed) in the CSR format of the library: a row-major one of another index type,
	// or a column-major one (transposed in parallel)
	template <typename Derived>
	CSRMatrix ToCSR(const Eigen::SparseCompressedBase<Derived>& M)
	{
		using Index = typename Derived::StorageIndex;
		static_assert(std::is_same<typename Derived::Scalar, double>::value, "fhhos4: matrices of doubles only");
		if (!M.isCompressed())
			throw Error("fhhos4::ToCSR(): the Eigen matrix is not compressed: call its makeCompressed() first.");
		if constexpr (Derived::IsRowMajor)
			return ToCSR(BasicCSRView<Index, Index>(M.rows(), M.cols(), M.outerIndexPtr(), M.innerIndexPtr(), M.valuePtr()));
		else // the arrays of the transpose
			return detail::Transpose(BasicCSRView<Index, Index>(M.cols(), M.rows(), M.outerIndexPtr(), M.innerIndexPtr(), M.valuePtr()));
	}

	// No copy: the CSR view of A^T, from the arrays of a column-major matrix A (compressed, with the index type
	// SparseIndex), which is A itself if A is symmetric, up to the rounding errors of its assembly (its asymmetry): for
	// the symmetric inputs of the setups (A_TT, A_FF, and A of a symmetric system) stored column-major. Not of a
	// temporary.
	template <typename Derived>
	CSRView TransposedView(const Eigen::SparseCompressedBase<Derived>& A)
	{
		static_assert(!Derived::IsRowMajor, "fhhos4::TransposedView(): for a column-major matrix (a row-major one converts to a CSRView itself)");
		static_assert(std::is_same<typename Derived::StorageIndex, SparseIndex>::value, "fhhos4::TransposedView(): the index type of the matrix is not fhhos4::SparseIndex: convert it with fhhos4::ToCSR() (a copy)");
		static_assert(std::is_same<typename Derived::Scalar, double>::value, "fhhos4: matrices of doubles only");
		if (!A.isCompressed())
			throw Error("fhhos4::TransposedView(): the Eigen matrix is not compressed: call its makeCompressed() first.");
		return CSRView(A.cols(), A.rows(), A.outerIndexPtr(), A.innerIndexPtr(), A.valuePtr());
	}

	template <typename Derived>
	CSRView TransposedView(const Eigen::SparseCompressedBase<Derived>&&) = delete;
}
