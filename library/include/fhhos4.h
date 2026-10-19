// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Pierre Matalon
#pragma once
#include <memory>
#include <type_traits>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "fhhos4_error.h"
#include "fhhos4_version.h"

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
// The class is compiled in the shared library fhhos4, with fhhos4's flags. Eigen objects cross its boundary: they must
// have the same memory layout in the calling code and in the library, which holds for Eigen 3.4 and 5.x, whatever the
// alignment settings (-march), with Eigen's default index type (checked by the constructor, which throws fhhos4::Error
// otherwise). The library never allocates memory that the calling code frees, nor the other way around.
namespace fhhos4
{
	// The sparse matrices of fhhos4: row-major, with the index type of its compilation (FHHOS4_SPARSE_INDEX: int by
	// default). The input matrices of this type, compressed, are read in place; the others are converted (copied).
	using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, FHHOS4_SPARSE_INDEX>;

	// The Krylov method preconditioned by one cycle of the AMG (Solver::Krylov)
	enum class Krylov
	{
		Default,  // FCG for a symmetric matrix (SetupSymmetric...()), GCR otherwise (SetupNonsymmetric...())
		FCG,      // flexible conjugate gradient (fhhos4 -s fcguamg), for a symmetric positive definite matrix
		GCR,      // GCR(10), generalized conjugate residual restarted every 10 iterations, Notay's economical version
		          // (fhhos4 -s gcruamg). Minimizes the residual norm: for any matrix
		BiCGSTAB, // Eigen's BiCGSTAB. Not flexible: it needs a linear preconditioner, a V- or W-cycle (not K)
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
	// smoothers are sequential, except the hybrid ones.
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
		                   // binding (OMP_PROC_BIND, OMP_PLACES) is set
		int Verbosity = 0; // 0: nothing; 1: one line per setup (number of levels, time) and per solve (iterations,
		                   // residual, time); 2: all the messages of fhhos4, with the residual at each iteration

		Solver()
			: Solver(EigenSample())
		{}
		~Solver();
		Solver(Solver&&) noexcept;
		Solver& operator=(Solver&&) noexcept;
		Solver(const Solver&) = delete;
		Solver& operator=(const Solver&) = delete;

		//----------- Uncondensed system: the interface of choice -----------//

		// Setup for a symmetric positive definite system, from its blocks (A_FT = A_TF^T), without the DoFs of the
		// Dirichlet faces. Computes the condensed matrix A = A_FF - A_FT A_TT^{-1} A_TF (CondensedMatrix()), and keeps what
		// Solve(b_T, b_F, x_T, x_F) needs: the blocks may be freed after the setup.
		// Checks the parameters, and the sizes and structure of the inputs against Dimension and the degrees.
		// - A_TT (cells x cells), A_TF (cells x faces), A_FF (faces x faces): the DoFs of each cell (resp. face) are
		//   contiguous, in the order of the rows (resp. columns) of A_TF. A_TT and A_FF are read from their lower
		//   triangular part (full or lower-only matrices).
		// - cellInterpOfOne, faceInterpOfOne: the interpolation of the function 1 on the polynomial bases of the cells and
		//   faces, one coefficient per cell and per face (in the order of the blocks of rows and columns of A_TF). The
		//   bases must be hierarchical (ordered by degree; not needed with PCoarsening::ToZero) with a constant first
		//   function phi_0: 1 = c phi_0, i.e. c = 1/phi_0. phi_0 may differ from one cell or face to another (e.g.
		//   L2-orthonormal bases: c = sqrt(|T|) on a cell T); c = 1 everywhere if phi_0 = 1.
		// Supported matrices: Eigen::SparseMatrix<double, Eigen::ColMajor or Eigen::RowMajor, int or long> (all of the
		// same type), compressed or not. Those of type SparseMatrix, compressed, are read in place; the others are
		// converted (copied).
		template <typename SparseMatrixType>
		void SetupSymmetric(const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FF,
		                    const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne);

		// Setup for a non-symmetric system (e.g. convection-diffusion), whose symmetric part is positive definite: as
		// SetupSymmetric(), with the block A_FT (faces x cells) and all the blocks read in full: the same
		// algorithm runs on the matrix itself. The couplings of the aggregation are computed from the averaged traces of
		// the blocks of A_TF and A_FT; the coarse operators are the Galerkin products P^T A P (as AGMG: Y. Notay, "An
		// aggregation-based algebraic multigrid method", ETNA 37, 2010); the Gauss-Seidel smoothers follow the downwind
		// order of each level (convection), in the pre- and post-smoothing; the K-cycle iterates with GCR. Krylov must not
		// be FCG, nor CoarseSolver Cholesky, nor the smoothers Chebyshev.
		// See README_library.md.
		template <typename SparseMatrixType>
		void SetupNonsymmetric(const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FT, const SparseMatrixType& A_FF,
		                       const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne);

		// Solves the system, after SetupSymmetric() or SetupNonsymmetric(): condenses the right-hand
		// side, b = b_F - A_FT A_TT^{-1} b_T, solves A x_F = b, then computes x_T = A_TT^{-1} (b_T - A_TF x_F).
		// x_F: the initial guess on input (e.g. zero), the solution on output; x_T: output. The result is that of
		// A x_F = b. Does not throw if the solver has not converged after MaxIterations iterations: see the result.
		Result Solve(const Eigen::Ref<const Eigen::VectorXd>& b_T, const Eigen::Ref<const Eigen::VectorXd>& b_F,
		             Eigen::Ref<Eigen::VectorXd> x_T, Eigen::Ref<Eigen::VectorXd> x_F);

		//---- Condensed system: to integrate the library into a code that condenses statically ----//

		// Setup for a symmetric positive definite matrix, from the condensed matrix and the blocks: as
		// SetupSymmetric(), with A (faces x faces), the condensed matrix A_FF - A_FT A_TT^{-1} A_TF, both
		// triangles, and A_FF optional (the algorithm only uses A_TT and A_TF, whose rows each belong to one cell).
		// The solver keeps a reference to A: A must stay alive and unchanged until the next setup or the destruction of
		// the solver (not a temporary). The blocks are only read by the setup: they may be freed after it.
		template <typename SparseMatrixType>
		void SetupSymmetricCondensed(const SparseMatrixType& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FF,
		                             const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne);

		template <typename SparseMatrixType>
		void SetupSymmetricCondensed(const SparseMatrixType& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF,
		                             const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne);

		// A temporary A would be destroyed before the solves
		template <typename SparseMatrixType>
		void SetupSymmetricCondensed(const SparseMatrixType&& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FF,
		                             const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne) = delete;
		template <typename SparseMatrixType>
		void SetupSymmetricCondensed(const SparseMatrixType&& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF,
		                             const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne) = delete;

		// Setup for a non-symmetric matrix: as SetupNonsymmetric(), with the condensed matrix A (see
		// SetupSymmetricCondensed()) and A_FF optional
		template <typename SparseMatrixType>
		void SetupNonsymmetricCondensed(const SparseMatrixType& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FT, const SparseMatrixType& A_FF,
		                                const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne);

		template <typename SparseMatrixType>
		void SetupNonsymmetricCondensed(const SparseMatrixType& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FT,
		                                const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne);

		// A temporary A would be destroyed before the solves
		template <typename SparseMatrixType>
		void SetupNonsymmetricCondensed(const SparseMatrixType&& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FT, const SparseMatrixType& A_FF,
		                                const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne) = delete;
		template <typename SparseMatrixType>
		void SetupNonsymmetricCondensed(const SparseMatrixType&& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FT,
		                                const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne) = delete;

		// Solves the condensed system A x = b, after any setup. x: the initial guess on input (e.g. zero), the solution on
		// output. Does not throw if the solver has not converged after MaxIterations iterations: see the result.
		Result Solve(const Eigen::Ref<const Eigen::VectorXd>& b, Eigen::Ref<Eigen::VectorXd> x);

		// Solves A x = b from x = 0
		Eigen::VectorXd Solve(const Eigen::Ref<const Eigen::VectorXd>& b)
		{
			Eigen::VectorXd x = Eigen::VectorXd::Zero(b.rows());
			Solve(b, x);
			return x;
		}

		//------------------------------------------------------------------//

		// The condensed matrix A of the last setup, both triangles: computed by SetupSymmetric() and SetupNonsymmetric(),
		// or given to SetupSymmetricCondensed() and SetupNonsymmetricCondensed() (converted to SparseMatrix if it has
		// another type). Valid until the next setup or the destruction of the solver.
		const SparseMatrix& CondensedMatrix() const;

		bool IsSetUp() const;
		int NumberOfLevels() const;

	private:
		// The Eigen objects of the interface (the sparse matrices, Eigen::Ref of the vectors) must have the same memory
		// layout in the calling code and in the library, compiled with their own versions of Eigen and flags. The
		// constructor checks it: it compares their sizes, then reads with the library's Eigen samples built by the
		// calling code.
		struct EigenLayout
		{
			int Version[3]; // of Eigen, for the error message
			int IndexSize, MatrixSize, ConstVectorSize, VectorSize; // sizeof(Eigen::Index), of SparseMatrix, of the Eigen::Ref
			bool VectorTriviallyCopyable; // Eigen::Ref<Eigen::VectorXd> is passed by value: how depends on it
		};

		static EigenLayout HostEigenLayout()
		{
			using ConstVector = Eigen::Ref<const Eigen::VectorXd>;
			using Vector = Eigen::Ref<Eigen::VectorXd>;
			return {
#ifdef EIGEN_PATCH_VERSION
				{ EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION, EIGEN_PATCH_VERSION }, // Eigen >= 5.0 (semantic versioning)
#else
				{ EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION },
#endif
				(int)sizeof(Eigen::Index), (int)sizeof(SparseMatrix), (int)sizeof(ConstVector), (int)sizeof(Vector),
				std::is_trivially_copyable<Vector>::value };
		}

		// The values of the samples: a 2 x 3 matrix with the non-zeros (0, 0), (0, 2), (1, 1), and a vector of size 3
		static double SampleValue(int i) { return 0.5 + i; }

		struct EigenSample
		{
			SparseMatrix Matrix;
			Eigen::VectorXd Vector;

			EigenSample()
				: Matrix(2, 3), Vector(3)
			{
				Matrix.insert(0, 0) = SampleValue(0);
				Matrix.insert(0, 2) = SampleValue(1);
				Matrix.insert(1, 1) = SampleValue(2);
				Matrix.makeCompressed();
				for (int i = 0; i < 3; i++)
					Vector[i] = SampleValue(i);
			}
		};

		Solver(EigenSample&& sample)
			: Solver(HostEigenLayout(), sample.Matrix, sample.Vector, sample.Vector)
		{}
		// sampleVector last: passed according to the calling code's layout, it is read once the layouts are checked
		Solver(const EigenLayout& host, const SparseMatrix& sampleMatrix, const Eigen::Ref<const Eigen::VectorXd>& sampleConstVector,
		       Eigen::Ref<Eigen::VectorXd> sampleVector);

		struct Impl;
		std::unique_ptr<Impl> _impl;
	};
}
