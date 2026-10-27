// The shared library fhhos4: the algebraic multigrid of fhhos4 for hybrid discretizations (U-AMG), for the codes
// that discretize the problem themselves (include/fhhos4.h). Compiled with fhhos4's flags and -fvisibility=hidden: only
// the classes of the public header are exported. The multigrid is built as in the program (SolverFactory,
// ApplyProgramArgumentDefaults), from the same code (AlgebraicSolverFactory, ApplyUncondensedAMGDefaults).
#include <chrono>
#include <limits>
#include <tuple>
#include <type_traits>
#include "fhhos4.h"
#include "ProgramArgumentsDefaults.h"
#include "Solver/AlgebraicSolverFactory.h"
#include "Solver/LibraryParameters.h"
#include "Solver/Krylov/PreconditionedBiCGSTAB.h"
#include "Solver/DistributedStaticCondensation.h"
#include "Solver/StaticCondensation.h"
#include "Utils/BorrowedSparseMatrix.h"
#include "Utils/SparseMatrixOps.h"

namespace fhhos4
{
	static_assert(is_same_v<fhhos4::SparseMatrix, ::SparseMatrix> && is_same_v<SparseIndex, SparseMatrixIndex>, "FHHOS4_SPARSE_INDEX (fhhos4_version.h) must be fhhos4's SparseMatrixIndex");

	using detail::BasicCSRView;
	using detail::CSRView;
	using detail::CSRMatrix;

	namespace
	{
		// Discards the messages (no stream buffer)
		ostream& NullStream()
		{
			static ostream stream(nullptr);
			return stream;
		}

		// During a call of the library: fhhos4's messages (Utils::Log()) according to the verbosity (on several processes,
		// those of the process of rank 0), the number of threads (Solver::Threads and Verbosity; the helpers: 0 and 0), and
		// the processes of the solver (Communicator::Current(): one by default). Restores them at the end, with the format
		// of cout, which printing the iterations changes.
		class CallScope
		{
		private:
			ostream* _savedLog;
			ios_base::fmtflags _savedCoutFlags;
			streamsize _savedCoutPrecision;
			int _savedThreads = 0;
			CurrentCommunicator _comm;
		public:
			// sharingProcesses: the processes of comm that run on the same CPUs as this one (Parallelism::SharingProcesses()),
			// among which the automatic number of threads divides the cores
			CallScope(const Solver& solver, const Communicator& comm, int sharingProcesses)
				: CallScope(solver.Threads, solver.Verbosity, comm, sharingProcesses)
			{}

			CallScope(int threads, int verbosity, const Communicator& comm = Communicator(), int sharingProcesses = 1)
				: _comm(comm)
			{
				_savedLog = Utils::LogStream;
				Utils::LogStream = verbosity >= 2 && comm.IsRoot() ? &cout : &NullStream();
				_savedCoutFlags = cout.flags();
				_savedCoutPrecision = cout.precision();
#ifdef _OPENMP
				if (threads <= 0)
					threads = Parallelism::WithoutHyperThreads(omp_get_max_threads(), sharingProcesses);
				if (threads != omp_get_max_threads())
				{
					_savedThreads = omp_get_max_threads();
					omp_set_num_threads(threads);
				}
#endif
			}

			~CallScope()
			{
				Utils::LogStream = _savedLog;
				cout.flags(_savedCoutFlags);
				cout.precision(_savedCoutPrecision);
#ifdef _OPENMP
				if (_savedThreads > 0)
					omp_set_num_threads(_savedThreads);
#endif
			}
		};

		double Seconds(chrono::steady_clock::time_point start)
		{
			return chrono::duration<double>(chrono::steady_clock::now() - start).count();
		}

		// The CSR format of a view: RowPointers[0] = 0 and non-decreasing, the columns in [0, Cols) and strictly increasing
		// in each row (no duplicates), the sizes and the number of non-zeros within the index type of the library. One
		// parallel pass on the indices.
		template <typename RowPointer, typename Column>
		void CheckCSR(const BasicCSRView<RowPointer, Column>& M, const string& name)
		{
			string matrix = "fhhos4: the matrix " + name + " (" + to_string(M.Rows) + " x " + to_string(M.Cols) + ")";
			if (M.Rows < 0 || M.Cols < 0)
				Utils::FatalError(matrix + " has a negative size.");
			if (!M.RowPointers)
			{
				if (M.Rows > 0)
					Utils::FatalError(matrix + " has no RowPointers.");
				return;
			}
			std::int64_t nnz = M.RowPointers[M.Rows];
			if (M.RowPointers[0] != 0 || nnz < 0)
				Utils::FatalError(matrix + " is not in the CSR format: RowPointers[0] = " + to_string(M.RowPointers[0]) + " (0 expected), RowPointers[Rows] = " + to_string(nnz) + ".");
			if (nnz > 0 && (!M.ColumnIndices || !M.Values))
				Utils::FatalError(matrix + " has non-zeros, but no ColumnIndices or Values.");
			constexpr std::int64_t maxIndex = numeric_limits<SparseMatrixIndex>::max();
			if (M.Rows > maxIndex || M.Cols > maxIndex || nnz > maxIndex)
				Utils::FatalError(matrix + " has too many rows, columns or non-zeros (" + to_string(nnz) + ") for the index type of the library (compile fhhos4 with -DSMALL_INDEX=OFF).");
			std::int64_t firstInvalidRow = M.Rows;
			#pragma omp parallel for reduction(min:firstInvalidRow) schedule(static)
			for (std::int64_t i = 0; i < M.Rows; ++i)
			{
				std::int64_t begin = M.RowPointers[i], end = M.RowPointers[i + 1];
				bool valid = begin <= end && end <= nnz;
				for (std::int64_t p = begin; valid && p < end; ++p)
				{
					std::int64_t col = M.ColumnIndices[p];
					valid = col >= 0 && col < M.Cols && (p == begin || col > (std::int64_t)M.ColumnIndices[p - 1]);
				}
				if (!valid)
					firstInvalidRow = min(firstInvalidRow, i);
			}
			if (firstInvalidRow < M.Rows)
				Utils::FatalError(matrix + " is not in the CSR format: in its row " + to_string(firstInvalidRow) + ", the columns are not strictly increasing (duplicates or unsorted: e.g. SciPy's sum_duplicates() sorts them and sums the duplicates), or outside [0, " + to_string(M.Cols) + "), or RowPointers decrease.");
		}

		// The matrix of a view, read in place, after CheckCSR()
		unique_ptr<BorrowedSparseMatrix> Borrow(const CSRView& M, const string& name)
		{
			CheckCSR(M, name);
			static const SparseMatrixIndex noRows[1] = { 0 }; // the RowPointers of a matrix without rows, if null
			return make_unique<BorrowedSparseMatrix>(M.Rows, M.Cols, M.RowPointers ? M.RowPointers : noRows, M.ColumnIndices, M.Values);
		}

		// Runs check() on all the processes of comm: if it throws on one of them, all of them throw (an error on one
		// process alone would leave the others waiting in the next collective call), with the message of the first
		// process where it threw
		template <class Check>
		void CollectiveCheck(const Communicator& comm, Check check)
		{
			if (!comm.IsDistributed())
			{
				check();
				return;
			}
			string message;
			try
			{
				check();
			}
			catch (const fhhos4::Error& e)
			{
				message = e.what();
			}
			int failedRank = comm.Min(message.empty() ? comm.Size() : comm.Rank());
			if (failedRank == comm.Size())
				return;
			vector<char> text(message.begin(), message.end());
			comm.Broadcast(text, failedRank);
			Utils::FatalError(string(text.begin(), text.end()) + (comm.Size() > 1 ? " (on the process " + to_string(failedRank) + ")" : ""));
		}
	}

	struct Solver::Impl
	{
		// Set by the setup
		bool Symmetric = true; // the matrix is symmetric (SetupSymmetric()), or not (SetupNonsymmetric())
		fhhos4::Krylov Krylov = fhhos4::Krylov::FCG;                       // Solver::Krylov, or its default for the symmetry of the matrix
		fhhos4::CoarseSolver CoarseSolver = fhhos4::CoarseSolver::Cholesky; // Solver::CoarseSolver, or its default for the symmetry of the matrix
		bool PolynomialSpaceQ = false; // Solver::PolynomialSpace
		int CellBlockSize = 0;
		int FaceBlockSize = 0;
		Eigen::Index NCellDoFs = 0;
		SparseMatrix A; // the condensed matrix, computed by the setup (on one process): the solvers keep a pointer to it
		bool Singular = false; // Solver::SingularMatrix, at the setup
		Vector Kernel;  // Singular: the interpolation of 1 on the faces (of the process), normalized (the kernel of A)
		unique_ptr<BorrowedSparseMatrix> GivenA_FT; // the caller's A_FT, read in place by Condensation (SetupNonsymmetric())
		unique_ptr<StaticCondensation> Condensation; // after a setup from the blocks: for Solve(b_T, b_F, x_T, x_F)
		unique_ptr<UncondensedAMG> Multigrid;
		unique_ptr<IterativeSolver> KrylovSolver; // preconditioned by Multigrid; null if Krylov is None
		bool IsSetUp = false;

		// The processes of the solver (Solver::Comm with MPI; one otherwise), and, on several processes, the faces of the
		// cells of the process (owned, then ghosts), the condensation and the condensed matrix, distributed
		Communicator Comm;
		int SharingProcesses = 1; // of Comm on the node of this process (see CallScope)
		shared_ptr<const IndexMap> Faces;
		unique_ptr<DistributedStaticCondensation> DistributedCondensation;
		DistributedMatrix DistributedA;
		Eigen::Index NFaceDoFs = 0; // of the process

		bool IsDistributed() const
		{
			return Comm.IsDistributed();
		}

		// The communicator of the solver, from Solver::Comm (mandatory with MPI)
		void SetCommunicator(const Solver& s)
		{
#if FHHOS4_MPI
			if (s.Comm == MPI_COMM_NULL)
				Utils::FatalError("fhhos4::Solver: set Comm (the library is compiled with MPI): MPI_COMM_SELF for the system of each process, or the communicator of the processes of a distributed system.");
			Comm = Communicator(s.Comm);
#else
			Comm = Communicator();
#endif
			SharingProcesses = Comm.SharingProcesses();
		}

		IterativeSolver* ActiveSolver() const
		{
			if (KrylovSolver)
				return KrylovSolver.get();
			return Multigrid.get();
		}

		// Checks the parameters, and computes the block sizes and the Krylov method for the symmetry of the matrix
		void CheckParameters(const Solver& s, bool symmetric)
		{
			if (s.Dimension != 2 && s.Dimension != 3)
				Utils::FatalError("fhhos4::Solver: set Dimension (2 or 3) before the setup, got " + to_string(s.Dimension) + ".");
			if (s.FaceDegree < 0)
				Utils::FatalError("fhhos4::Solver: set FaceDegree (>= 0) before the setup, got " + to_string(s.FaceDegree) + ".");
			if (s.CellDegree < 0)
				Utils::FatalError("fhhos4::Solver: set CellDegree (>= 0) before the setup, got " + to_string(s.CellDegree) + ".");
			switch (s.Krylov)
			{
				case fhhos4::Krylov::Default:
					Krylov = symmetric ? fhhos4::Krylov::FCG : fhhos4::Krylov::GCR;
					break;
				case fhhos4::Krylov::FCG:
				case fhhos4::Krylov::GCR:
				case fhhos4::Krylov::BiCGSTAB:
				case fhhos4::Krylov::None:
					Krylov = s.Krylov;
					break;
				default:
					Utils::FatalError("fhhos4::Solver: invalid Krylov (" + to_string((int)s.Krylov) + ").");
			}
			if (Krylov == fhhos4::Krylov::BiCGSTAB && s.Cycle == fhhos4::Cycle::K)
				Utils::FatalError("fhhos4::Solver: BiCGSTAB is not flexible: it needs a linear preconditioner, a V- or W-cycle, not the K-cycle (non-linear: it solves the coarse levels with FCG or GCR iterations). Set Cycle to Cycle::V or Cycle::W, or Krylov to Krylov::FCG or Krylov::GCR.");
			Symmetric = symmetric;
			CoarseSolver = s.CoarseSolver != fhhos4::CoarseSolver::Default ? s.CoarseSolver : (symmetric ? fhhos4::CoarseSolver::Cholesky : fhhos4::CoarseSolver::LU);
			if (!symmetric && Krylov == fhhos4::Krylov::FCG)
				Utils::FatalError("fhhos4::Solver: the flexible conjugate gradient requires a symmetric matrix: set Krylov to Krylov::GCR (or Default) for a non-symmetric one.");
			if (!symmetric && CoarseSolver == fhhos4::CoarseSolver::Cholesky)
				Utils::FatalError("fhhos4::Solver: the Cholesky factorization requires a symmetric matrix: set CoarseSolver to CoarseSolver::LU (or Default) for a non-symmetric one.");
			if (!symmetric && s.SingularMatrix)
				Utils::FatalError("fhhos4::Solver: SingularMatrix requires a symmetric matrix (SetupSymmetric()): the kernel of the transpose of a non-symmetric matrix is not the interpolation of 1 in general.");
			if (s.PreSmoothingIterations < -1 || s.PostSmoothingIterations < -1)
				Utils::FatalError("fhhos4::Solver: PreSmoothingIterations and PostSmoothingIterations must be >= 0, or -1 (the smoother's default), got " + to_string(s.PreSmoothingIterations) + " and " + to_string(s.PostSmoothingIterations) + ".");
			for (fhhos4::Smoother smoother : { s.PreSmoother, s.PostSmoother })
			{
				if (!symmetric && SmootherFactory::IsChebyshev(LibraryParameters::ProgramSmoother(smoother)))
					Utils::FatalError("fhhos4::Solver: the Chebyshev smoothers require a symmetric positive definite matrix: leave PreSmoother and PostSmoother Default for a non-symmetric one.");
			}
			PolynomialSpaceQ = LibraryParameters::ProgramPolynomialSpace(s.PolynomialSpace);
			CellBlockSize = Utils::PolynomialSpaceDimension(s.Dimension, s.CellDegree, PolynomialSpaceQ);
			FaceBlockSize = Utils::PolynomialSpaceDimension(s.Dimension - 1, s.FaceDegree, PolynomialSpaceQ);
		}

		// The sizes and structure of the inputs, against the dimension and the degrees. A_FT: null if the matrix is
		// symmetric.
		void CheckInputs(const Solver& s, const SparseMatrix& A_TT, const SparseMatrix& A_TF, const SparseMatrix* A_FT, const SparseMatrix& A_FF, const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne) const
		{
			string space = string(", PolynomialSpace ") + (PolynomialSpaceQ ? "Q" : "P");
			string cellDoFs = to_string(CellBlockSize) + " DoFs per cell (CellDegree " + to_string(s.CellDegree) + ", Dimension " + to_string(s.Dimension) + space + ")";
			string faceDoFs = to_string(FaceBlockSize) + " DoFs per face (FaceDegree " + to_string(s.FaceDegree) + ", Dimension " + to_string(s.Dimension) + space + ")";
			Eigen::Index nCellDoFs = A_TF.rows();
			Eigen::Index nFaceDoFs = A_TF.cols();
			auto size = [](const SparseMatrix& M) { return to_string(M.rows()) + " x " + to_string(M.cols()); };
			if (nCellDoFs % CellBlockSize != 0)
				Utils::FatalError("fhhos4::Solver: the number of rows of A_TF (" + to_string(nCellDoFs) + ") is not a multiple of the " + cellDoFs + ".");
			if (nFaceDoFs % FaceBlockSize != 0)
				Utils::FatalError("fhhos4::Solver: the number of columns of A_TF (" + to_string(nFaceDoFs) + ") is not a multiple of the " + faceDoFs + ".");
			if (A_TT.rows() != nCellDoFs || A_TT.cols() != nCellDoFs)
				Utils::FatalError("fhhos4::Solver: A_TT is " + size(A_TT) + ", expected " + to_string(nCellDoFs) + " x " + to_string(nCellDoFs) + " (the rows of A_TF).");
			if (A_FT && (A_FT->rows() != nFaceDoFs || A_FT->cols() != nCellDoFs))
				Utils::FatalError("fhhos4::Solver: A_FT is " + size(*A_FT) + ", expected " + to_string(nFaceDoFs) + " x " + to_string(nCellDoFs) + " (the transpose of A_TF).");
			if (A_FF.rows() != nFaceDoFs || A_FF.cols() != nFaceDoFs)
				Utils::FatalError("fhhos4::Solver: A_FF is " + size(A_FF) + ", expected " + to_string(nFaceDoFs) + " x " + to_string(nFaceDoFs) + " (the columns of A_TF).");
			if (cellInterpOfOne.rows() != nCellDoFs / CellBlockSize || faceInterpOfOne.rows() != nFaceDoFs / FaceBlockSize)
				Utils::FatalError("fhhos4::Solver: the interpolations of 1 on the cells and faces must have one coefficient per cell (" + to_string(nCellDoFs / CellBlockSize) + ": the rows of A_TF, with " + cellDoFs + ") and per face (" + to_string(nFaceDoFs / FaceBlockSize) + ": the columns of A_TF, with " + faceDoFs + "), got " + to_string(cellInterpOfOne.rows()) + " and " + to_string(faceInterpOfOne.rows()) + ".");

			// A_TT is block diagonal
			for (BigNumber i = 0; i < A_TT.rows(); i++)
			{
				for (SparseMatrix::InnerIterator it(A_TT, i); it; ++it)
				{
					if (it.col() / CellBlockSize != i / CellBlockSize && it.value() != 0)
						Utils::FatalError("fhhos4::Solver: A_TT has a non-zero (" + to_string(i) + ", " + to_string(it.col()) + ") outside the diagonal blocks of the cells, with " + cellDoFs + ".");
				}
			}
			CheckInterpOfOne(cellInterpOfOne, "cell");
			CheckInterpOfOne(faceInterpOfOne, "face");
		}

		static void CheckInterpOfOne(const Eigen::Ref<const Eigen::VectorXd>& interpOfOne, const string& entity)
		{
			for (Eigen::Index i = 0; i < interpOfOne.rows(); i++)
			{
				if (!std::isfinite(interpOfOne[i]) || interpOfOne[i] == 0)
					Utils::FatalError("fhhos4::Solver: the interpolation of 1 on the " + entity + " bases is not finite and non-zero on the " + entity + " " + to_string(i) + ".");
			}
		}

		// The pre-smoother when Solver::PreSmoother is empty, whatever the Krylov method (PERFORMANCE.md): for a symmetric
		// matrix, the block Chebyshev of the fourth kind, whatever the number of threads (the user's choice: the same
		// results on any number of threads; on 1 thread, slower than the block Gauss-Seidel); for a non-symmetric one
		// (Chebyshev requires a symmetric matrix), the hybrid block Gauss-Seidel on several threads, the block Gauss-Seidel
		// on 1
		SmootherType DefaultPreSmoother() const
		{
			if (Symmetric)
				return SmootherType::BlockChebyshevFourthKind;
			return Parallelism::NThreads() > 1 ? SmootherType::HybridBlockGaussSeidel : SmootherType::BlockGaussSeidel;
		}

		// The number of iterations of the smoother when Solver::Pre/PostSmoothingIterations is -1: the degree 2 for the
		// Chebyshev smoothers (PERFORMANCE.md), 1 otherwise
		static int DefaultSmoothingIterations(SmootherType smoother)
		{
			return SmootherFactory::IsChebyshev(smoother) ? 2 : 1;
		}

		// The arguments of the program fhhos4 that give the same multigrid
		ProgramArguments Arguments(const Solver& s) const
		{
			ProgramArguments args;
			args.Solver.Tolerance = s.Tolerance;
			args.Solver.MaxIterations = s.MaxIterations;
			args.Solver.PrintIterationResults = s.Verbosity >= 2;
			MultigridArguments& mg = args.Solver.MG;
			mg.Cycle = LibraryParameters::ProgramCycle(s.Cycle);
			mg.WLoops = s.Cycle == fhhos4::Cycle::W ? 2 : 1;
			mg.P_CS = LibraryParameters::ProgramPCoarsening(s.PCoarsening);
			mg.CoarseningProlongationCode = (int)LibraryParameters::ProgramProlongation(s.Prolongation);
			mg.KCycleKrylov = Krylov == fhhos4::Krylov::GCR || !Symmetric ? KrylovMethod::GCR : KrylovMethod::FCG; // GCR under GCR (Notay)
			mg.PreSmoother = s.PreSmoother != fhhos4::Smoother::Default ? LibraryParameters::ProgramSmoother(s.PreSmoother) : DefaultPreSmoother();
			mg.PostSmoother = LibraryParameters::ProgramSmoother(s.PostSmoother); // Default: the program's, the adjoint of the pre-smoother (rbgs after bgs, the same Chebyshev...)
			mg.PreSmoothingIterations = s.PreSmoothingIterations >= 0 ? s.PreSmoothingIterations : DefaultSmoothingIterations(mg.PreSmoother);
			// PostSmoother Default: of the type of the pre-smoother
			mg.PostSmoothingIterations = s.PostSmoothingIterations >= 0 ? s.PostSmoothingIterations : DefaultSmoothingIterations(mg.PostSmoother != SmootherType::Default ? mg.PostSmoother : mg.PreSmoother);
			mg.CoarseningFactor = s.CoarseningFactor;
			mg.MatrixMaxSizeForCoarsestLevel = s.CoarseMatrixMaxSize;
			mg.CoarseSolverCode = LibraryParameters::ProgramCoarseSolver(CoarseSolver, Symmetric);
			mg.UAMGMixedPrecision = s.MixedPrecision;
			mg.UAMGBlockCSRStorage = s.BlockCSRStorage;
			bool defaultCoarseOperator = true, defaultCycle = false, defaultHPCoarseningStgy = true;
			ApplyUncondensedAMGDefaults(mg, s.FaceDegree, defaultCoarseOperator, defaultCycle, defaultHPCoarseningStgy);
			return args;
		}

		// A_FT: null if the matrix is symmetric; kept until the next setup (the solves read it)
		void Setup(const Solver& s, const SparseMatrix& A_TT, const SparseMatrix& A_TF, unique_ptr<BorrowedSparseMatrix> A_FT, const SparseMatrix& A_FF, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne)
		{
			auto start = chrono::steady_clock::now();
			Reset();
			GivenA_FT = std::move(A_FT);
			NCellDoFs = A_TF.rows();
			NFaceDoFs = A_TF.cols();

			ProgramArguments args = Arguments(s);
			CreateSolvers(args, s);

			Condensation = make_unique<StaticCondensation>(A_TT, A_TF, GivenA_FT.get(), CellBlockSize);
			A = Condensation->CondensedMatrix(A_FF, Symmetric);
			if (s.SingularMatrix)
				SetKernel(faceInterpOfOne);

			// The blocks are only read during the setup (UncondensedAMG::Setup())
			SparseMatrix symmetric; // empty A_FT: symmetric matrix
			ActiveSolver()->Setup(A, A_TT, A_TF, GivenA_FT ? *GivenA_FT : symmetric, A_FF, cellInterpOfOne, faceInterpOfOne);
			IsSetUp = true;

			if (s.Verbosity == 1)
				cout << "fhhos4: setup of " << Multigrid->NumberOfLevels() << " levels in " << Seconds(start) << " s" << endl;
		}

		// On several processes: from the blocks of the uncondensed system, the rows of the process (see
		// Solver::SetupSymmetric()), A_FF in full; A_FT on the rows of all the faces of faces and the columns of the cells
		// of the process (DistributedStaticCondensation), empty if the system is symmetric
		void SetupDistributed(const Solver& s, const SparseMatrix& A_TT, const SparseMatrix& A_TF, SparseMatrix A_FT, const SparseMatrix& A_FF, shared_ptr<const IndexMap> faces, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne)
		{
			auto start = chrono::steady_clock::now();
			Reset();
			Faces = faces;
			NCellDoFs = A_TF.rows();
			NFaceDoFs = A_FF.rows();

			ProgramArguments args = Arguments(s);
			CreateSolvers(args, s);

			DistributedCondensation = make_unique<DistributedStaticCondensation>(A_TT, A_TF, std::move(A_FT), faces, CellBlockSize, FaceBlockSize);
			DistributedA = DistributedCondensation->CondensedMatrix(A_FF);
			if (s.SingularMatrix)
				SetKernel(faceInterpOfOne);
			ActiveSolver()->Setup(DistributedA, A_TT, A_TF, DistributedCondensation->A_FT(), faces, cellInterpOfOne, faceInterpOfOne);
			IsSetUp = true;

			if (s.Verbosity == 1 && Comm.IsRoot())
				cout << "fhhos4: setup of " << Multigrid->NumberOfLevels() << " levels on " << Comm.Size() << " processes in " << Seconds(start) << " s" << endl;
		}

		// Solver::SingularMatrix: the kernel of A, the interpolation of 1 on the faces (c on the first DoF of each face, 0
		// on the others; on several processes, those of the process), checked: A z = 0 up to the rounding errors.
		// Otherwise the projections of the right-hand sides would silently solve another problem (e.g. a Dirichlet face in
		// the system, or a wrong interpolation of 1).
		void SetKernel(const Vector& faceInterpOfOne)
		{
			Singular = true;
			Kernel = SingularSolver::Kernel(faceInterpOfOne, FaceBlockSize, NFaceDoFs);
			double relative;
			if (IsDistributed())
				relative = Parallel::Norm(DistributedA * Kernel) / sqrt(Comm.Sum(DistributedA.Diag.squaredNorm() + DistributedA.Offd.squaredNorm()));
			else
				relative = Vector(A.selfadjointView<Eigen::Lower>() * Kernel).norm() / A.norm(); // A: lower triangular part
			if (!(relative < 1e-8))
				Utils::FatalError("fhhos4::Solver: SingularMatrix, but the interpolation of 1 is not in the kernel of the matrix (|A z| / (|A| |z|) = " + to_string(relative) + ", z: the interpolation of 1 on the faces): is there a Dirichlet face in the system, or is the interpolation of 1 wrong?");
		}

		void Reset()
		{
			IsSetUp = false;
			Singular = false;
			Kernel.resize(0);
			KrylovSolver.reset();
			Multigrid.reset();
			Condensation.reset();
			DistributedCondensation.reset();
			DistributedA = DistributedMatrix();
			Faces.reset();
			A = SparseMatrix();
		}

		// The multigrid, and the Krylov method preconditioned by one cycle (SolverPreconditioner), as SolverFactory (-s
		// fcguamg, gcruamg), except the residual replacement: the tolerance applies to the true residual
		// (Result::RelativeResidual = |b - Ax|/|b|)
		void CreateSolvers(const ProgramArguments& args, const Solver& s)
		{
			ExportModule out;
			Multigrid.reset(AlgebraicSolverFactory::CreateUncondensedAMG(args, s.Dimension, s.FaceDegree, PolynomialSpaceQ, CellBlockSize, FaceBlockSize, FaceBlockSize, out, AlgebraicSolverFactory::AlgebraicCoarseSolvers(s.Dimension, out)));
			Multigrid->SingularMatrix = s.SingularMatrix;
			if (Krylov == fhhos4::Krylov::FCG)
			{
				auto fcg = make_unique<FlexibleConjugateGradient>(1);
				fcg->ResidualRecomputation = 20;
				fcg->Precond = SolverPreconditioner(Multigrid.get());
				KrylovSolver = std::move(fcg);
			}
			else if (Krylov == fhhos4::Krylov::GCR)
			{
				auto gcr = make_unique<GeneralizedConjugateResidual>();
				gcr->ResidualReplacement = true;
				gcr->Precond = SolverPreconditioner(Multigrid.get());
				KrylovSolver = std::move(gcr);
			}
			else if (Krylov == fhhos4::Krylov::BiCGSTAB)
			{
				auto bicgstab = make_unique<PreconditionedBiCGSTAB>();
				bicgstab->Precond = SolverPreconditioner(Multigrid.get());
				KrylovSolver = std::move(bicgstab);
			}
			ActiveSolver()->PrintIterationResults = args.Solver.PrintIterationResults;
		}

		void CheckSetUp(const char* function) const
		{
			if (!IsSetUp)
				Utils::FatalError(string("fhhos4::Solver::") + function + "(): a setup function must be called first (and succeed).");
		}

		void CheckVectorSize(Eigen::Index size, const char* vector) const
		{
			if (size != NFaceDoFs)
				Utils::FatalError(string("fhhos4::Solver: the vector ") + vector + " has " + to_string(size) + " rows, expected " + to_string(NFaceDoFs) + " (the faces DoFs" + (IsDistributed() ? " of the process" : "") + ").");
		}

		void CheckCellVectorSize(Eigen::Index size, const char* vector) const
		{
			if (size != NCellDoFs)
				Utils::FatalError(string("fhhos4::Solver: the vector ") + vector + " has " + to_string(size) + " rows, expected " + to_string(NCellDoFs) + " (the cells DoFs" + (IsDistributed() ? " of the process" : "") + ").");
		}
	};

	Solver::Solver()
		: _impl(make_unique<Impl>())
	{}

	Solver::~Solver() = default;
	Solver::Solver(Solver&&) noexcept = default;
	Solver& Solver::operator=(Solver&&) noexcept = default;

	void Solver::Setup(bool symmetric, const CSRView& A_TT, const CSRView& A_TF, const CSRView* A_FT, const CSRView& A_FF,
	                   detail::ConstVectorView cellInterpOfOne, detail::ConstVectorView faceInterpOfOne, const GhostNumbers& faceGhosts,
	                   const GhostNumbers& cellGhosts)
	{
		_impl->SetCommunicator(*this);
		CallScope scope(*this, _impl->Comm, _impl->SharingProcesses);
		if (_impl->IsDistributed())
		{
			SetupDistributed(symmetric, A_TT, A_TF, A_FT, A_FF, cellInterpOfOne, faceInterpOfOne, faceGhosts, cellGhosts);
			return;
		}
		if (faceGhosts.Count != 0 || cellGhosts.Count != 0)
			Utils::FatalError("fhhos4::Solver: faceGhosts or cellGhosts on one process (Comm): the ghosts are the faces and cells of other processes.");
		_impl->CheckParameters(*this, symmetric);
		// The inputs, read in place, checked before the previous setup is replaced
		unique_ptr<BorrowedSparseMatrix> a_TT = Borrow(A_TT, "A_TT"), a_TF = Borrow(A_TF, "A_TF");
		unique_ptr<BorrowedSparseMatrix> a_FT = A_FT ? Borrow(*A_FT, "A_FT") : nullptr;
		unique_ptr<BorrowedSparseMatrix> a_FF = Borrow(A_FF, "A_FF");
		Eigen::Map<const Vector> cellOne(cellInterpOfOne.Data, cellInterpOfOne.Size), faceOne(faceInterpOfOne.Data, faceInterpOfOne.Size);
		_impl->CheckInputs(*this, *a_TT, *a_TF, a_FT.get(), *a_FF, cellOne, faceOne);
		_impl->Setup(*this, *a_TT, *a_TF, std::move(a_FT), *a_FF, cellOne, faceOne);
	}

	// On several processes, the parameters are those of the process of rank 0 (the same on all of them)
	void Solver::SetupDistributed(bool symmetric, const CSRView& A_TT, const CSRView& A_TF, const CSRView* A_FT, const CSRView& A_FF,
	                              detail::ConstVectorView cellInterpOfOne, detail::ConstVectorView faceInterpOfOne, const GhostNumbers& faceGhosts,
	                              const GhostNumbers& cellGhosts)
	{
		const Communicator& comm = _impl->Comm;
		unique_ptr<BorrowedSparseMatrix> a_TT, a_TF, a_FT, a_FF;
		shared_ptr<const IndexMap> faces;
		CollectiveCheck(comm, [&]()
		{
			_impl->CheckParameters(*this, symmetric);
			for (auto [ghosts, name] : { pair(&faceGhosts, "faceGhosts"), pair(&cellGhosts, "cellGhosts") })
				if (ghosts->Count < 0 || (ghosts->Count > 0 && !ghosts->Numbers))
					Utils::FatalError(string("fhhos4::Solver: ") + name + ": a negative Count, or no Numbers.");
			if (symmetric && cellGhosts.Count != 0)
				Utils::FatalError("fhhos4::Solver: cellGhosts are the ghost cells of A_FT, for a non-symmetric system only (SetupNonsymmetric()).");
			a_TT = Borrow(A_TT, "A_TT");
			a_TF = Borrow(A_TF, "A_TF");
			if (A_FT)
				a_FT = Borrow(*A_FT, "A_FT");
			a_FF = Borrow(A_FF, "A_FF");
		});
		int cbs = _impl->CellBlockSize, fbs = _impl->FaceBlockSize;
		CollectiveCheck(comm, [&]()
		{
			auto size = [](const SparseMatrix& M) { return to_string(M.rows()) + " x " + to_string(M.cols()); };
			if (a_TF->rows() % cbs != 0 || a_TT->rows() != a_TF->rows() || a_TT->cols() != a_TT->rows())
				Utils::FatalError("fhhos4::Solver: on several processes, A_TT (" + size(*a_TT) + ") and A_TF (" + size(*a_TF) + ") have the rows of the cells of the process (" + to_string(cbs) + " DoFs per cell), A_TT the same columns.");
			if (a_FF->rows() % fbs != 0 || a_FF->cols() != a_TF->cols() || a_FF->cols() != a_FF->rows() + faceGhosts.Count * fbs)
				Utils::FatalError("fhhos4::Solver: on several processes, A_FF (" + size(*a_FF) + ") has the rows of the faces of the process (" + to_string(fbs) + " DoFs per face), and A_FF and A_TF (" + size(*a_TF) + ") the columns of the faces of the process then of its " + to_string(faceGhosts.Count) + " ghosts (faceGhosts).");
			if (a_FT && (a_FT->rows() != a_FF->rows() || a_FT->cols() != a_TF->rows() + cellGhosts.Count * cbs))
				Utils::FatalError("fhhos4::Solver: on several processes, A_FT (" + size(*a_FT) + ") has the rows of the faces of the process (those of A_FF, " + size(*a_FF) + "), and the columns of the cells of the process (the rows of A_TF, " + size(*a_TF) + ") then of its " + to_string(cellGhosts.Count) + " ghost cells (cellGhosts).");
			if (cellInterpOfOne.Size != a_TF->rows() / cbs || faceInterpOfOne.Size != a_FF->rows() / fbs)
				Utils::FatalError("fhhos4::Solver: on several processes, the interpolations of 1 have the cells (" + to_string(a_TF->rows() / cbs) + ") and faces (" + to_string(a_FF->rows() / fbs) + ") of the process, got " + to_string(cellInterpOfOne.Size) + " and " + to_string(faceInterpOfOne.Size) + ".");
			// A_TT is block diagonal
			for (BigNumber i = 0; i < a_TT->rows(); i++)
				for (SparseMatrix::InnerIterator it(*a_TT, i); it; ++it)
					if (it.col() / cbs != i / cbs && it.value() != 0)
						Utils::FatalError("fhhos4::Solver: A_TT has a non-zero (" + to_string(i) + ", " + to_string(it.col()) + ") outside the diagonal blocks of the cells.");
			Impl::CheckInterpOfOne(Eigen::Map<const Vector>(cellInterpOfOne.Data, cellInterpOfOne.Size), "cell");
			Impl::CheckInterpOfOne(Eigen::Map<const Vector>(faceInterpOfOne.Data, faceInterpOfOne.Size), "face");
		});
		// The faces of the process and its ghosts (collective, checks the ghosts on all the processes)
		IndexMap owners(comm, a_FF->rows() / fbs);
		faces = make_shared<const IndexMap>(IndexMap::InOrder(owners, vector<BigNumber>(faceGhosts.Numbers, faceGhosts.Numbers + faceGhosts.Count)));
		// Non-symmetric: A_FT on the rows of all the faces of the process's cells (those of other processes fetched from
		// them) and the columns of these cells, what the condensation of these cells and their aggregation read (collective,
		// checks the ghost cells on all the processes)
		SparseMatrix localA_FT;
		if (!symmetric)
		{
			IndexMap cellOwners(comm, a_TF->rows() / cbs);
			auto cells = make_shared<const IndexMap>(IndexMap::InOrder(cellOwners, vector<BigNumber>(cellGhosts.Numbers, cellGhosts.Numbers + cellGhosts.Count)));
			localA_FT = DistributedMatrix::FromLocalColumns(make_shared<const IndexMap>(owners), cells, *a_FT, fbs, cbs).RowsWithGhosts(*faces, cellOwners);
		}
		Eigen::Map<const Vector> cellOne(cellInterpOfOne.Data, cellInterpOfOne.Size), faceOne(faceInterpOfOne.Data, faceInterpOfOne.Size);
		_impl->SetupDistributed(*this, *a_TT, *a_TF, std::move(localA_FT), *a_FF, faces, cellOne, faceOne);
	}

	Result Solver::SolveCondensed(detail::ConstVectorView b, detail::VectorView x)
	{
		_impl->CheckSetUp("Solve");
		CallScope scope(*this, _impl->Comm, _impl->SharingProcesses);
		CollectiveCheck(_impl->Comm, [&]()
		{
			_impl->CheckVectorSize(b.Size, "b");
			_impl->CheckVectorSize(x.Size, "x");
		});
		auto start = chrono::steady_clock::now();

		Result result;
		Vector rhs = Eigen::Map<const Vector>(b.Data, b.Size);
		if (_impl->Singular)
		{
			// The right-hand side projected onto the image of A (orthogonal to its kernel)
			double component = Parallel::Dot(rhs, _impl->Kernel), norm = Parallel::Norm(rhs);
			if (Verbosity >= 1 && _impl->Comm.IsRoot() && abs(component) > 1e-10 * norm)
				cout << "fhhos4: the right-hand side is not in the image of the singular matrix: its component along the kernel removed (relative size " << abs(component) / norm << ")" << endl;
			rhs -= component * _impl->Kernel;
		}
		Eigen::Map<Vector> callerX(x.Data, x.Size);
		if (Parallel::Norm(rhs) == 0)
		{
			callerX.setZero();
			result.Converged = true;
			return result;
		}

		Vector solution = callerX;
		bool xEquals0 = _impl->Comm.All((solution.array() == 0).all()); // skips the products by the initial guess, as the program does
		IterativeSolver* solver = _impl->ActiveSolver();
		solver->Tolerance = Tolerance;
		solver->MaxIterations = MaxIterations;
		solver->Solve(rhs, solution, xEquals0);
		callerX = solution;

		result.Iterations = solver->IterationCount;
		result.RelativeResidual = solver->LastIterationResult.NormalizedResidualNorm;
		result.Converged = result.RelativeResidual < Tolerance;
		if (Verbosity == 1 && _impl->Comm.IsRoot())
			cout << "fhhos4: " << result.Iterations << " iterations, relative residual " << result.RelativeResidual << (result.Converged ? "" : " (not converged)") << ", " << Seconds(start) << " s" << endl;
		return result;
	}

	Result Solver::SolveUncondensed(detail::ConstVectorView b_T, detail::ConstVectorView b_F, detail::VectorView x_T, detail::VectorView x_F)
	{
		_impl->CheckSetUp("Solve");
		CallScope scope(*this, _impl->Comm, _impl->SharingProcesses);
		CollectiveCheck(_impl->Comm, [&]()
		{
			_impl->CheckCellVectorSize(b_T.Size, "b_T");
			_impl->CheckCellVectorSize(x_T.Size, "x_T");
			_impl->CheckVectorSize(b_F.Size, "b_F");
			_impl->CheckVectorSize(x_F.Size, "x_F");
		});

		if (_impl->DistributedCondensation)
		{
			const DistributedStaticCondensation& condensation = *_impl->DistributedCondensation;
			Vector invA_TT_b_T = condensation.SolveCells(Eigen::Map<const Vector>(b_T.Data, b_T.Size));
			Vector b = condensation.CondensedRHS(Eigen::Map<const Vector>(b_F.Data, b_F.Size), invA_TT_b_T);
			Result result = SolveCondensed({ b.data(), b.size() }, x_F);
			Eigen::Map<Vector>(x_T.Data, x_T.Size) = condensation.CellUnknowns(invA_TT_b_T, Eigen::Map<const Vector>(x_F.Data, x_F.Size));
			return result;
		}

		const StaticCondensation& condensation = *_impl->Condensation;
		Vector invA_TT_b_T = condensation.SolveCells(Eigen::Map<const Vector>(b_T.Data, b_T.Size));
		Vector b = condensation.CondensedRHS(Eigen::Map<const Vector>(b_F.Data, b_F.Size), invA_TT_b_T);
		Result result = SolveCondensed({ b.data(), b.size() }, x_F);
		Eigen::Map<Vector>(x_T.Data, x_T.Size) = condensation.CellUnknowns(invA_TT_b_T, Eigen::Map<const Vector>(x_F.Data, x_F.Size));
		return result;
	}

	CSRView Solver::CondensedMatrixView() const
	{
		_impl->CheckSetUp("CondensedMatrix");
		if (_impl->IsDistributed())
			Utils::FatalError("fhhos4::Solver::CondensedMatrix(): not available on several processes.");
		const SparseMatrix& A = _impl->A;
		return CSRView(A.rows(), A.cols(), A.outerIndexPtr(), A.innerIndexPtr(), A.valuePtr());
	}

	bool Solver::IsSetUp() const
	{
		return _impl->IsSetUp;
	}

	int Solver::NumberOfLevels() const
	{
		_impl->CheckSetUp("NumberOfLevels");
		return _impl->Multigrid->NumberOfLevels();
	}

	//------------------------------------------------------------------------------------------------------------------//
	//                       Helpers: the blocks from the global matrix of the system (copies)                        //
	//------------------------------------------------------------------------------------------------------------------//

	namespace detail
	{
		struct CSRMatrix::Impl
		{
			SparseMatrix Matrix;
		};

		CSRMatrix::CSRMatrix() = default;
		CSRMatrix::CSRMatrix(unique_ptr<Impl> impl)
			: _impl(std::move(impl))
		{}
		CSRMatrix::~CSRMatrix() = default;
		CSRMatrix::CSRMatrix(CSRMatrix&&) noexcept = default;
		CSRMatrix& CSRMatrix::operator=(CSRMatrix&&) noexcept = default;

		CSRView CSRMatrix::View() const&
		{
			if (!_impl)
				return CSRView();
			const SparseMatrix& M = _impl->Matrix;
			return CSRView(M.rows(), M.cols(), M.outerIndexPtr(), M.innerIndexPtr(), M.valuePtr());
		}
	}

	namespace
	{
		CSRMatrix Own(SparseMatrix&& M)
		{
			auto impl = make_unique<CSRMatrix::Impl>();
			impl->Matrix = std::move(M);
			impl->Matrix.makeCompressed();
			return CSRMatrix(std::move(impl));
		}

		// The lists: Count >= 0, Offsets[0] = 0 and non-decreasing, the DoFs in [0, n). Returns their number of DoFs.
		template <typename Index>
		BigNumber CheckDoFLists(const DoFLists<Index>& lists, std::int64_t n, const string& name)
		{
			string error = "fhhos4: the DoF lists of the " + name + " (" + to_string(lists.Count) + ")";
			if (lists.Count < 0)
				Utils::FatalError(error + ": negative Count.");
			if (lists.Count == 0)
				return 0;
			if (!lists.Offsets)
				Utils::FatalError(error + " have no Offsets.");
			if (lists.Offsets[0] != 0)
				Utils::FatalError(error + ": Offsets[0] = " + to_string(lists.Offsets[0]) + " (0 expected).");
			for (std::int64_t i = 0; i < lists.Count; ++i)
			{
				if (lists.Offsets[i + 1] < lists.Offsets[i])
					Utils::FatalError(error + ": the Offsets decrease at " + to_string(i) + ".");
			}
			std::int64_t nDoFs = lists.Offsets[lists.Count];
			if (nDoFs > 0 && !lists.DoFs)
				Utils::FatalError(error + " have no DoFs.");
			if (nDoFs > numeric_limits<SparseMatrixIndex>::max())
				Utils::FatalError(error + " have too many DoFs (" + to_string(nDoFs) + ") for the index type of the library (compile fhhos4 with -DSMALL_INDEX=OFF).");
			std::int64_t firstInvalid = nDoFs;
			#pragma omp parallel for reduction(min:firstInvalid) schedule(static)
			for (std::int64_t p = 0; p < nDoFs; ++p)
			{
				if (lists.DoFs[p] < 0 || lists.DoFs[p] >= n)
					firstInvalid = min(firstInvalid, p);
			}
			if (firstInvalid < nDoFs)
				Utils::FatalError(error + ": the DoF " + to_string(lists.DoFs[firstInvalid]) + " (DoFs[" + to_string(firstInvalid) + "]) is outside [0, " + to_string(n) + ").");
			return nDoFs;
		}

		// The DoFs of the global numbering: of a cell, of a face, or of none
		enum class DoFKind : uint8_t { None, Cell, Face };

		// The rows of M whose global indices are rows[0..nRows), in this order, restricted to the columns of the kind
		// columnKind, renumbered by their positions in their lists. lower: only the columns whose position is <= that of
		// the row (the lower triangular part of a diagonal block). M: after CheckCSR().
		template <typename RowPointer, typename Column, typename Index>
		SparseMatrix ExtractBlock(const BasicCSRView<RowPointer, Column>& M, const Index* rows, BigNumber nRows, const vector<SparseMatrixIndex>& position, const vector<DoFKind>& kind, DoFKind columnKind, BigNumber nCols, bool lower)
		{
			ThreadLocal<vector<pair<SparseMatrixIndex, double>>> rowsOfThreads;
			return SparseMatrixOps::BuildByRows(nRows, nCols, 1, [&](BigNumber firstRow, BigNumber lastRow, SparseMatrixOps::RowChunk& chunk)
			{
				vector<pair<SparseMatrixIndex, double>>& row = rowsOfThreads.Local();
				for (BigNumber i = firstRow; i < lastRow; ++i)
				{
					std::int64_t g = rows[i];
					row.clear();
					for (std::int64_t p = M.RowPointers[g]; p < (std::int64_t)M.RowPointers[g + 1]; ++p)
					{
						std::int64_t col = M.ColumnIndices[p];
						if (kind[col] == columnKind && (!lower || position[col] <= (SparseMatrixIndex)i))
							row.push_back({ position[col], M.Values[p] });
					}
					auto byColumn = [](const pair<SparseMatrixIndex, double>& x, const pair<SparseMatrixIndex, double>& y) { return x.first < y.first; };
					if (!is_sorted(row.begin(), row.end(), byColumn)) // the lists don't follow the global numbering
						sort(row.begin(), row.end(), byColumn);
					for (const pair<SparseMatrixIndex, double>& coeff : row)
						chunk.Add(coeff.first, coeff.second);
					chunk.EndRow();
				}
			});
		}

		template <typename RowPointer, typename Column, typename Index>
		detail::HybridBlocks ExtractBlocks(const BasicCSRView<RowPointer, Column>& M, const DoFLists<Index>& cells, const DoFLists<Index>& faces, bool symmetric)
		{
			CallScope scope(0, 0);
			if (M.Rows != M.Cols)
				Utils::FatalError("fhhos4: the global matrix M (" + to_string(M.Rows) + " x " + to_string(M.Cols) + ") must be square.");
			CheckCSR(M, "M");
			BigNumber nT = CheckDoFLists(cells, M.Rows, "cells");
			BigNumber nF = CheckDoFLists(faces, M.Rows, "faces");

			// The position of each DoF of M in the cells (resp. faces) DoFs: in the order of the lists
			vector<SparseMatrixIndex> position(M.Rows, -1);
			vector<DoFKind> kind(M.Rows, DoFKind::None);
			for (auto [lists, nDoFs, k] : { tuple{ &cells, nT, DoFKind::Cell }, { &faces, nF, DoFKind::Face } })
			{
				for (BigNumber p = 0; p < nDoFs; ++p)
				{
					std::int64_t g = lists->DoFs[p];
					if (kind[g] != DoFKind::None)
						Utils::FatalError("fhhos4: the DoF " + to_string(g) + " of the global matrix is in two lists (of the cells or faces).");
					kind[g] = k;
					position[g] = p;
				}
			}

			detail::HybridBlocks blocks;
			blocks.A_TT = Own(ExtractBlock(M, cells.DoFs, nT, position, kind, DoFKind::Cell, nT, symmetric));
			blocks.A_TF = Own(ExtractBlock(M, cells.DoFs, nT, position, kind, DoFKind::Face, nF, false));
			if (!symmetric)
				blocks.A_FT = Own(ExtractBlock(M, faces.DoFs, nF, position, kind, DoFKind::Cell, nT, false));
			blocks.A_FF = Own(ExtractBlock(M, faces.DoFs, nF, position, kind, DoFKind::Face, nF, symmetric));
			return blocks;
		}

		void CheckVectorSize(std::int64_t size, std::int64_t expected, const string& vector, const string& of)
		{
			if (size != expected)
				Utils::FatalError("fhhos4: the vector " + vector + " has " + to_string(size) + " rows, expected " + to_string(expected) + " (" + of + ").");
		}
	}

	namespace detail
	{
		template <typename RowPointer, typename Column, typename Index>
		HybridBlocks ExtractSymmetricBlocks(const BasicCSRView<RowPointer, Column>& M, const DoFLists<Index>& cells, const DoFLists<Index>& faces)
		{
			return ExtractBlocks(M, cells, faces, true);
		}

		template <typename RowPointer, typename Column, typename Index>
		HybridBlocks ExtractNonsymmetricBlocks(const BasicCSRView<RowPointer, Column>& M, const DoFLists<Index>& cells, const DoFLists<Index>& faces)
		{
			return ExtractBlocks(M, cells, faces, false);
		}

		template <typename RowPointer, typename Column>
		CSRMatrix ToCSR(const BasicCSRView<RowPointer, Column>& M)
		{
			CallScope scope(0, 0);
			CheckCSR(M, "M");
			SparseMatrix copy(M.Rows, M.Cols);
			BigNumber nnz = M.NonZeros();
			copy.resizeNonZeros(nnz);
			SparseMatrixIndex* rowPointers = copy.outerIndexPtr();
			SparseMatrixIndex* columns = copy.innerIndexPtr();
			double* values = copy.valuePtr();
			#pragma omp parallel for schedule(static)
			for (BigNumber i = 0; i <= (BigNumber)M.Rows; ++i)
				rowPointers[i] = M.RowPointers[i];
			#pragma omp parallel for schedule(static)
			for (BigNumber p = 0; p < nnz; ++p)
			{
				columns[p] = M.ColumnIndices[p];
				values[p] = M.Values[p];
			}
			return Own(std::move(copy));
		}

		template <typename Index>
		CSRMatrix Transpose(const BasicCSRView<Index, Index>& M)
		{
			CallScope scope(0, 0);
			CheckCSR(M, "M");
			Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor, Index>> arrays(M.Rows, M.Cols, M.NonZeros(), M.RowPointers, M.ColumnIndices, M.Values);
			return Own(SparseMatrixOps::Transpose(arrays));
		}

		template <typename Index>
		void Gather(ConstVectorView x, const DoFLists<Index>& cells, const DoFLists<Index>& faces, VectorView x_T, VectorView x_F)
		{
			CallScope scope(0, 0);
			BigNumber nT = CheckDoFLists(cells, x.Size, "cells");
			BigNumber nF = CheckDoFLists(faces, x.Size, "faces");
			CheckVectorSize(x_T.Size, nT, "x_T", "the DoFs of the cells");
			CheckVectorSize(x_F.Size, nF, "x_F", "the DoFs of the faces");
			#pragma omp parallel for schedule(static)
			for (BigNumber p = 0; p < nT; ++p)
				x_T.Data[p] = x.Data[cells.DoFs[p]];
			#pragma omp parallel for schedule(static)
			for (BigNumber p = 0; p < nF; ++p)
				x_F.Data[p] = x.Data[faces.DoFs[p]];
		}

		template <typename Index>
		void Scatter(ConstVectorView x_T, ConstVectorView x_F, const DoFLists<Index>& cells, const DoFLists<Index>& faces, VectorView x)
		{
			CallScope scope(0, 0);
			BigNumber nT = CheckDoFLists(cells, x.Size, "cells");
			BigNumber nF = CheckDoFLists(faces, x.Size, "faces");
			CheckVectorSize(x_T.Size, nT, "x_T", "the DoFs of the cells");
			CheckVectorSize(x_F.Size, nF, "x_F", "the DoFs of the faces");
			#pragma omp parallel for schedule(static)
			for (BigNumber p = 0; p < nT; ++p)
				x.Data[cells.DoFs[p]] = x_T.Data[p];
			#pragma omp parallel for schedule(static)
			for (BigNumber p = 0; p < nF; ++p)
				x.Data[faces.DoFs[p]] = x_F.Data[p];
		}
	}

	// The index types of the helpers: std::int32_t or std::int64_t for each array
#define FHHOS4_INSTANTIATE_EXTRACTION(RowPointer, Column, Index) \
	template detail::HybridBlocks detail::ExtractSymmetricBlocks<RowPointer, Column, Index>(const BasicCSRView<RowPointer, Column>&, const DoFLists<Index>&, const DoFLists<Index>&); \
	template detail::HybridBlocks detail::ExtractNonsymmetricBlocks<RowPointer, Column, Index>(const BasicCSRView<RowPointer, Column>&, const DoFLists<Index>&, const DoFLists<Index>&);
#define FHHOS4_INSTANTIATE_CSR(RowPointer, Column) \
	template CSRMatrix detail::ToCSR<RowPointer, Column>(const BasicCSRView<RowPointer, Column>&); \
	FHHOS4_INSTANTIATE_EXTRACTION(RowPointer, Column, std::int32_t) \
	FHHOS4_INSTANTIATE_EXTRACTION(RowPointer, Column, std::int64_t)

	FHHOS4_INSTANTIATE_CSR(std::int32_t, std::int32_t)
	FHHOS4_INSTANTIATE_CSR(std::int32_t, std::int64_t)
	FHHOS4_INSTANTIATE_CSR(std::int64_t, std::int32_t)
	FHHOS4_INSTANTIATE_CSR(std::int64_t, std::int64_t)

	template CSRMatrix detail::Transpose<std::int32_t>(const BasicCSRView<std::int32_t, std::int32_t>&);
	template CSRMatrix detail::Transpose<std::int64_t>(const BasicCSRView<std::int64_t, std::int64_t>&);
	template void detail::Gather<std::int32_t>(detail::ConstVectorView, const DoFLists<std::int32_t>&, const DoFLists<std::int32_t>&, detail::VectorView, detail::VectorView);
	template void detail::Gather<std::int64_t>(detail::ConstVectorView, const DoFLists<std::int64_t>&, const DoFLists<std::int64_t>&, detail::VectorView, detail::VectorView);
	template void detail::Scatter<std::int32_t>(detail::ConstVectorView, detail::ConstVectorView, const DoFLists<std::int32_t>&, const DoFLists<std::int32_t>&, detail::VectorView);
	template void detail::Scatter<std::int64_t>(detail::ConstVectorView, detail::ConstVectorView, const DoFLists<std::int64_t>&, const DoFLists<std::int64_t>&, detail::VectorView);
}
