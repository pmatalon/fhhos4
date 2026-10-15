// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 Pierre Matalon
#pragma once
#include <memory>
#include <string>
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
// The class is compiled in the shared library fhhos4, with fhhos4's flags. Eigen objects cross
// its boundary: the calling code must use the same version of Eigen and the same alignment settings (checked by the
// constructor, which throws fhhos4::Error otherwise). The library never allocates memory that the calling code frees,
// nor the other way around.
namespace fhhos4
{
	// The sparse matrices of fhhos4: row-major, with the index type of its compilation (FHHOS4_SPARSE_INDEX: int by
	// default). The input matrices of this type, compressed, are read in place; the others are converted (copied).
	using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, FHHOS4_SPARSE_INDEX>;

	struct Result
	{
		int Iterations = 0;
		double RelativeResidual = 0; // |b - A x| / |b|, in the condensed system
		bool Converged = false;      // RelativeResidual < Tolerance
	};

	// Set the parameters, then call a setup function, then Solve(). The parameters are read by the setup (Tolerance and
	// MaxIterations by Solve()). The defaults are those of the program fhhos4 (-s fcguamg), used in the papers, except the
	// smoothers on several threads (PreSmoother, and the number of its iterations).
	// Not thread-safe: one object must not be used by several threads at the same time.
	class FHHOS4_API Solver
	{
	public:
		//------------- Discretization: to set -------------//

		int Dimension = -1;  // 2 or 3
		int FaceDegree = -1; // k: degree of the polynomials on the faces
		int CellDegree = -1; // l: degree of the polynomials in the cells (k - 1, k or k + 1 in HHO)

		//--------------------- Solver ---------------------//

		// The Krylov method preconditioned by one cycle of the AMG:
		// - "" (default): "fcg" for a symmetric matrix (SetupSymmetric...()), "gcr" otherwise (SetupNonsymmetric...());
		// - "fcg": flexible conjugate gradient (fhhos4 -s fcguamg), for a symmetric positive definite matrix;
		// - "gcr": GCR(10), generalized conjugate residual restarted every 10 iterations, Notay's economical version
		//   (fhhos4 -s gcruamg). Minimizes the residual norm: for any matrix;
		// - "bicgstab": Eigen's BiCGSTAB. Not flexible: it needs a linear preconditioner, a V- or W-cycle (not K);
		// - "none": the cycles of the AMG alone, as a stationary iterative method (fhhos4 -s uamg): slower at k >= 1.
		std::string Krylov;
		double Tolerance = 1e-8; // stopping criterion of Solve(): |b - A x| / |b| < Tolerance
		int MaxIterations = 200; // Solve() stops after MaxIterations iterations, converged or not

		//-------------------- Multigrid -------------------//

		char Cycle = 'K';                  // 'V', 'W' or 'K' (K-cycle: each coarse level is solved by 2 iterations of FCG,
		                                   // or of GCR with Krylov "gcr" or a non-symmetric matrix)
		// Number of iterations of the pre- and post-smoothers (>= 0, 0: none); -1 (default): the smoother's, 2 for the
		// Chebyshev smoothers (their degree), 1 for the others
		int PreSmoothingIterations = -1;
		int PostSmoothingIterations = -1;
		// Smoothers (codes of fhhos4 -smoothers, see fhhos4 -h). PreSmoother empty (default): chosen by the setup:
		// - for a symmetric matrix (SetupSymmetric...()), "bcheb4": the block Chebyshev of the fourth kind, a polynomial in
		//   D^-1 A (D: the block diagonal, one block per face) whose degree is the number of smoothing iterations (each
		//   one, a block Jacobi iteration with its own damping factor): parallel, and its results don't depend on the
		//   number of threads ("bcheb": the first kind; symmetric matrices only);
		// - for a non-symmetric matrix, on several threads, "hbgs": the hybrid block Gauss-Seidel (Gauss-Seidel in the
		//   rows of each thread, Jacobi between the threads: parallel, except on the levels where the rows of the threads
		//   are poorly separated; the results then depend on the number of threads, not on the run); on 1 thread, "bgs":
		//   the block Gauss-Seidel (sequential: "hbgs" on 1 thread), one block per face.
		// PostSmoother empty (default): the backward version of PreSmoother ("rbgs" after "bgs", "hrbgs" after "hbgs"),
		// or PreSmoother itself if it is symmetric (Chebyshev).
		// Non-symmetric matrix (SetupNonsymmetric...()): on a level with one-sided couplings (convection), the Gauss-Seidel
		// smoothers keep their type, block size and number of sweeps, but not their direction: all their sweeps
		// (backward and symmetric ones included) follow the downwind order of the level. The other smoothers (Jacobi,
		// SOR with a relaxation parameter) ignore it.
		std::string PreSmoother;
		std::string PostSmoother;
		double CoarseningFactor = 3.8;     // target ratio of the numbers of cells of two successive h-levels
		int CoarseMatrixMaxSize = 1000;    // the coarsening stops when the matrix has at most this many rows
		std::string CoarseSolver;          // solver of the coarsest level (fhhos4 -coarse-solver); "" (default): sparse
		                                   // Cholesky ("ch") for a symmetric matrix, sparse LU ("lu") otherwise

		//-------------------- Execution -------------------//

		int Threads = 0;   // number of OpenMP threads during the calls of the library; 0: the current OpenMP setting, but
		                   // at most one thread per physical core (no hyper-threading) unless OMP_NUM_THREADS or OpenMP's
		                   // binding (OMP_PROC_BIND, OMP_PLACES) is set
		int Verbosity = 0; // 0: nothing; 1: one line per setup (number of levels, time) and per solve (iterations,
		                   // residual, time); 2: all the messages of fhhos4, with the residual at each iteration

		Solver()
			: Solver(HostEigenConfiguration())
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
		//   bases must be hierarchical (ordered by degree) with a constant first function phi_0: 1 = c phi_0, i.e.
		//   c = 1/phi_0. phi_0 may differ from one cell or face to another (e.g. L2-orthonormal bases: c = sqrt(|T|) on a
		//   cell T); c = 1 everywhere if phi_0 = 1.
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
		// be "fcg", nor CoarseSolver "ch".
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
		// The configuration of Eigen that determines the memory layout of its objects
		struct EigenConfiguration
		{
			int World, Major, Minor, MaxAlignBytes, IndexSize;
		};

		static EigenConfiguration HostEigenConfiguration()
		{
			return { EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION, EIGEN_MAX_ALIGN_BYTES, (int)sizeof(Eigen::Index) };
		}

		Solver(const EigenConfiguration& host);

		struct Impl;
		std::unique_ptr<Impl> _impl;
	};
}
