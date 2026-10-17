// The shared library fhhos4: the algebraic multigrid of fhhos4 for hybrid discretizations (U-AMG), for the codes
// that discretize the problem themselves (include/fhhos4.h). Compiled with fhhos4's flags and -fvisibility=hidden: only
// the classes of the public header are exported. The multigrid is built as in the program (SolverFactory,
// ApplyProgramArgumentDefaults), from the same code (AlgebraicSolverFactory, ApplyUncondensedAMGDefaults).
#include <chrono>
#include <limits>
#include <type_traits>
#include "fhhos4.h"
#include "ProgramArgumentsDefaults.h"
#include "Solver/AlgebraicSolverFactory.h"
#include "Solver/LibraryParameters.h"
#include "Solver/Krylov/PreconditionedBiCGSTAB.h"
#include "Solver/StaticCondensation.h"
#include "Utils/SparseMatrixOps.h"

namespace fhhos4
{
	static_assert(is_same_v<fhhos4::SparseMatrix, ::SparseMatrix>, "FHHOS4_SPARSE_INDEX (fhhos4_version.h) must be fhhos4's SparseMatrixIndex");

	namespace
	{
		// Discards the messages (no stream buffer)
		ostream& NullStream()
		{
			static ostream stream(nullptr);
			return stream;
		}

		// During a call of the library: fhhos4's messages (Utils::Log()) according to the verbosity, and the number of
		// threads of the parameters. Restores them at the end, with the format of cout, which printing the iterations changes.
		class CallScope
		{
		private:
			ostream* _savedLog;
			ios_base::fmtflags _savedCoutFlags;
			streamsize _savedCoutPrecision;
			int _savedThreads = 0;
		public:
			CallScope(const Solver& solver)
			{
				_savedLog = Utils::LogStream;
				Utils::LogStream = solver.Verbosity >= 2 ? &cout : &NullStream();
				_savedCoutFlags = cout.flags();
				_savedCoutPrecision = cout.precision();
#ifdef _OPENMP
				int threads = solver.Threads > 0 ? solver.Threads : Parallelism::WithoutHyperThreads(omp_get_max_threads());
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

		// A matrix of the calling code in fhhos4's format (row-major, fhhos4's index type, compressed): the caller's matrix
		// itself if it has this format (no copy), a converted copy otherwise. Not copyable: it may point to its copy.
		class InputMatrix
		{
		private:
			SparseMatrix _copy; // empty if the caller's matrix is used
			const SparseMatrix* _matrix;
		public:
			template <typename SparseMatrixType>
			InputMatrix(const SparseMatrixType& M, const string& name)
			{
				if constexpr (is_same_v<SparseMatrixType, SparseMatrix>)
				{
					if (M.isCompressed())
					{
						_matrix = &M;
						return;
					}
				}
				if (M.nonZeros() > (Eigen::Index)numeric_limits<SparseMatrixIndex>::max())
					Utils::FatalError("fhhos4::Solver: " + name + " has too many non-zeros for the index type of fhhos4 (compile fhhos4 with -DSMALL_INDEX=OFF).");
				if constexpr (!SparseMatrixType::IsRowMajor && is_same_v<typename SparseMatrixType::StorageIndex, SparseMatrixIndex>)
				{
					if (M.isCompressed())
					{
						// The arrays of a compressed column-major matrix are those of its transpose in row-major, transposed
						// again in parallel: the matrix of Eigen's conversion, which is sequential
						_copy = SparseMatrixOps::Transpose(Eigen::Map<const SparseMatrix>(M.cols(), M.rows(), M.nonZeros(), M.outerIndexPtr(), M.innerIndexPtr(), M.valuePtr()));
						_matrix = &_copy;
						return;
					}
				}
				_copy = M;
				_copy.makeCompressed();
				_matrix = &_copy;
			}
			InputMatrix(const InputMatrix&) = delete;
			InputMatrix& operator=(const InputMatrix&) = delete;

			const SparseMatrix& Get() const
			{
				return *_matrix;
			}

			// The matrix to keep after the setup: the caller's matrix, or the copy, moved into storage
			const SparseMatrix* Keep(SparseMatrix& storage)
			{
				if (_matrix == &_copy)
				{
					storage = std::move(_copy);
					_matrix = &storage;
				}
				return _matrix;
			}
		};
	}

	struct Solver::Impl
	{
		// Set by the setup
		bool Symmetric = true; // the matrix is symmetric (SetupSymmetric...()), or not (SetupNonsymmetric...())
		fhhos4::Krylov Krylov = fhhos4::Krylov::FCG;                       // Solver::Krylov, or its default for the symmetry of the matrix
		fhhos4::CoarseSolver CoarseSolver = fhhos4::CoarseSolver::Cholesky; // Solver::CoarseSolver, or its default for the symmetry of the matrix
		int CellBlockSize = 0;
		int FaceBlockSize = 0;
		Eigen::Index NCellDoFs = 0;
		const SparseMatrix* A = nullptr; // the solvers keep a pointer to it: the caller's matrix, or OwnA
		SparseMatrix OwnA;               // A if the library computes or converts it, empty otherwise
		unique_ptr<StaticCondensation> Condensation; // after a setup from the blocks: for Solve(b_T, b_F, x_T, x_F)
		unique_ptr<UncondensedAMG> Multigrid;
		unique_ptr<IterativeSolver> KrylovSolver; // preconditioned by Multigrid; null if Krylov is None
		bool IsSetUp = false;

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
			if (s.PreSmoothingIterations < -1 || s.PostSmoothingIterations < -1)
				Utils::FatalError("fhhos4::Solver: PreSmoothingIterations and PostSmoothingIterations must be >= 0, or -1 (the smoother's default), got " + to_string(s.PreSmoothingIterations) + " and " + to_string(s.PostSmoothingIterations) + ".");
			for (fhhos4::Smoother smoother : { s.PreSmoother, s.PostSmoother })
			{
				if (!symmetric && SmootherFactory::IsChebyshev(LibraryParameters::ProgramSmoother(smoother)))
					Utils::FatalError("fhhos4::Solver: the Chebyshev smoothers require a symmetric positive definite matrix: leave PreSmoother and PostSmoother Default for a non-symmetric one.");
			}
			CellBlockSize = Utils::Binomial(s.CellDegree + s.Dimension, s.CellDegree);
			FaceBlockSize = Utils::Binomial(s.FaceDegree + s.Dimension - 1, s.FaceDegree);
		}

		// The sizes and structure of the inputs, against the dimension and the degrees. A, A_FF: null if not given; A_FT:
		// null if the matrix is symmetric.
		void CheckInputs(const Solver& s, const SparseMatrix* A, const SparseMatrix& A_TT, const SparseMatrix& A_TF, const SparseMatrix* A_FT, const SparseMatrix* A_FF, const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne) const
		{
			string cellDoFs = to_string(CellBlockSize) + " DoFs per cell (CellDegree " + to_string(s.CellDegree) + ", Dimension " + to_string(s.Dimension) + ")";
			string faceDoFs = to_string(FaceBlockSize) + " DoFs per face (FaceDegree " + to_string(s.FaceDegree) + ", Dimension " + to_string(s.Dimension) + ")";
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
			if (A_FF && (A_FF->rows() != nFaceDoFs || A_FF->cols() != nFaceDoFs))
				Utils::FatalError("fhhos4::Solver: A_FF is " + size(*A_FF) + ", expected " + to_string(nFaceDoFs) + " x " + to_string(nFaceDoFs) + " (the columns of A_TF).");
			if (A && (A->rows() != nFaceDoFs || A->cols() != nFaceDoFs))
				Utils::FatalError("fhhos4::Solver: A is " + size(*A) + ", expected " + to_string(nFaceDoFs) + " x " + to_string(nFaceDoFs) + " (the columns of A_TF).");
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
			mg.KCycleKrylov = Krylov == fhhos4::Krylov::GCR || !Symmetric ? KrylovMethod::GCR : KrylovMethod::FCG; // GCR under GCR (Notay)
			mg.PreSmoother = s.PreSmoother != fhhos4::Smoother::Default ? LibraryParameters::ProgramSmoother(s.PreSmoother) : DefaultPreSmoother();
			mg.PostSmoother = LibraryParameters::ProgramSmoother(s.PostSmoother); // Default: the program's, the adjoint of the pre-smoother (rbgs after bgs, the same Chebyshev...)
			mg.PreSmoothingIterations = s.PreSmoothingIterations >= 0 ? s.PreSmoothingIterations : DefaultSmoothingIterations(mg.PreSmoother);
			// PostSmoother Default: of the type of the pre-smoother
			mg.PostSmoothingIterations = s.PostSmoothingIterations >= 0 ? s.PostSmoothingIterations : DefaultSmoothingIterations(mg.PostSmoother != SmootherType::Default ? mg.PostSmoother : mg.PreSmoother);
			mg.CoarseningFactor = s.CoarseningFactor;
			mg.MatrixMaxSizeForCoarsestLevel = s.CoarseMatrixMaxSize;
			mg.CoarseSolverCode = LibraryParameters::ProgramCoarseSolver(CoarseSolver);
			bool defaultCoarseOperator = true, defaultCycle = false, defaultHPCoarseningStgy = true;
			ApplyUncondensedAMGDefaults(mg, s.FaceDegree, defaultCoarseOperator, defaultCycle, defaultHPCoarseningStgy);
			return args;
		}

		// A: kept until the next setup; null to compute it from the blocks. A_FT: null if the matrix is symmetric. A_FF: null
		// if not given (not with A null).
		void Setup(const Solver& s, InputMatrix* A, const SparseMatrix& A_TT, const SparseMatrix& A_TF, const SparseMatrix* A_FT, const SparseMatrix* A_FF, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne)
		{
			auto start = chrono::steady_clock::now();
			IsSetUp = false;
			KrylovSolver.reset();
			Multigrid.reset();
			Condensation.reset();
			this->A = nullptr;
			OwnA = SparseMatrix();
			NCellDoFs = A_TF.rows();

			ProgramArguments args = Arguments(s);
			ExportModule out;
			Multigrid.reset(AlgebraicSolverFactory::CreateUncondensedAMG(args, s.Dimension, s.FaceDegree, CellBlockSize, FaceBlockSize, FaceBlockSize, out, AlgebraicSolverFactory::AlgebraicCoarseSolvers(s.Dimension, out)));

			if (A)
				this->A = A->Keep(OwnA);
			else
			{
				Condensation = make_unique<StaticCondensation>(A_TT, A_TF, A_FT, CellBlockSize);
				OwnA = Condensation->CondensedMatrix(*A_FF, Symmetric);
				this->A = &OwnA;
			}
			SparseMatrix computedA_FF; // empty if A_FF is given, or not needed (the algorithm of the paper only uses A_TT and A_TF)
			if (!A_FF && Multigrid->A_F_FNeeded())
				computedA_FF = StaticCondensation(A_TT, A_TF, A_FT, CellBlockSize).FaceBlock(*this->A);

			// The Krylov method preconditioned by one cycle (SolverPreconditioner), as SolverFactory (-s fcguamg, gcruamg),
			// except the residual replacement: the tolerance applies to the true residual (Result::RelativeResidual =
			// |b - Ax|/|b|)
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

			// The blocks are only read during the setup (UncondensedAMG::Setup())
			SparseMatrix symmetric; // empty A_FT: symmetric matrix
			ActiveSolver()->Setup(*this->A, A_TT, A_TF, A_FT ? *A_FT : symmetric, A_FF ? *A_FF : computedA_FF, cellInterpOfOne, faceInterpOfOne);
			IsSetUp = true;

			if (s.Verbosity == 1)
				cout << "fhhos4: setup of " << Multigrid->NumberOfLevels() << " levels in " << Seconds(start) << " s" << endl;
		}

		void CheckSetUp(const char* function) const
		{
			if (!IsSetUp)
				Utils::FatalError(string("fhhos4::Solver::") + function + "(): a setup function must be called first (and succeed).");
		}

		void CheckVectorSize(Eigen::Index size, const char* vector) const
		{
			if (size != A->rows())
				Utils::FatalError(string("fhhos4::Solver: the vector ") + vector + " has " + to_string(size) + " rows, expected " + to_string(A->rows()) + " (the faces DoFs).");
		}

		void CheckCellVectorSize(Eigen::Index size, const char* vector) const
		{
			if (size != NCellDoFs)
				Utils::FatalError(string("fhhos4::Solver: the vector ") + vector + " has " + to_string(size) + " rows, expected " + to_string(NCellDoFs) + " (the cells DoFs).");
		}
	};

	Solver::Solver(const EigenLayout& host, const SparseMatrix& sampleMatrix, const Eigen::Ref<const Eigen::VectorXd>& sampleConstVector,
	               Eigen::Ref<Eigen::VectorXd> sampleVector)
	{
		EigenLayout library = HostEigenLayout(); // the layout of the library's compilation
		auto version = [](const EigenLayout& l) { return to_string(l.Version[0]) + "." + to_string(l.Version[1]) + "." + to_string(l.Version[2]); };
		auto sizes = [](const EigenLayout& l)
		{
			return to_string(l.IndexSize) + ", " + to_string(l.MatrixSize) + ", " + to_string(l.ConstVectorSize) + ", " + to_string(l.VectorSize) + (l.VectorTriviallyCopyable ? ", trivially copyable" : "");
		};
		string incompatible = "fhhos4::Solver: the Eigen objects of the interface have different memory layouts in the calling code (Eigen " + version(host) + ") and in the library fhhos4 (Eigen " + version(library) + "): ";
		string advice = ". Use Eigen 3.4 or 5.x, with its default index type (EIGEN_DEFAULT_DENSE_INDEX_TYPE).";
		if (host.IndexSize != library.IndexSize || host.MatrixSize != library.MatrixSize || host.ConstVectorSize != library.ConstVectorSize ||
		    host.VectorSize != library.VectorSize || host.VectorTriviallyCopyable != library.VectorTriviallyCopyable)
			throw Error(incompatible + "sizes of Eigen::Index, SparseMatrix, Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>: " + sizes(host) + " in the calling code, " + sizes(library) + " in the library" + advice);

		// The same sizes: the samples (fhhos4.h), whose dimensions are checked before their arrays are read
		const auto* outer = sampleMatrix.outerIndexPtr();
		const auto* inner = sampleMatrix.innerIndexPtr();
		const double* values = sampleMatrix.valuePtr();
		bool sameMatrix = sampleMatrix.rows() == 2 && sampleMatrix.cols() == 3 && sampleMatrix.isCompressed() && sampleMatrix.nonZeros() == 3 &&
		                  outer[0] == 0 && outer[1] == 2 && outer[2] == 3 && inner[0] == 0 && inner[1] == 2 && inner[2] == 1 &&
		                  values[0] == SampleValue(0) && values[1] == SampleValue(1) && values[2] == SampleValue(2);
		bool sameVectors = sampleConstVector.size() == 3 && sampleVector.size() == 3 && sampleConstVector.data() == sampleVector.data();
		for (int i = 0; sameVectors && i < 3; i++)
			sameVectors = sampleConstVector[i] == SampleValue(i) && sampleVector[i] == SampleValue(i);
		if (!sameMatrix || !sameVectors)
			throw Error(incompatible + "the library reads a different " + (sameMatrix ? "vector" : "sparse matrix") + " than the calling code built" + advice);

		_impl = make_unique<Impl>();
	}

	Solver::~Solver() = default;
	Solver::Solver(Solver&&) noexcept = default;
	Solver& Solver::operator=(Solver&&) noexcept = default;

	template <typename SparseMatrixType>
	void Solver::SetupSymmetricCondensed(const SparseMatrixType& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FF,
	                                     const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
	{
		CallScope scope(*this);
		_impl->CheckParameters(*this, true);
		InputMatrix a(A, "A"), a_TT(A_TT, "A_TT"), a_TF(A_TF, "A_TF"), a_FF(A_FF, "A_FF");
		_impl->CheckInputs(*this, &a.Get(), a_TT.Get(), a_TF.Get(), nullptr, &a_FF.Get(), cellInterpOfOne, faceInterpOfOne);
		_impl->Setup(*this, &a, a_TT.Get(), a_TF.Get(), nullptr, &a_FF.Get(), cellInterpOfOne, faceInterpOfOne);
	}

	template <typename SparseMatrixType>
	void Solver::SetupSymmetricCondensed(const SparseMatrixType& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF,
	                                     const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
	{
		CallScope scope(*this);
		_impl->CheckParameters(*this, true);
		InputMatrix a(A, "A"), a_TT(A_TT, "A_TT"), a_TF(A_TF, "A_TF");
		_impl->CheckInputs(*this, &a.Get(), a_TT.Get(), a_TF.Get(), nullptr, nullptr, cellInterpOfOne, faceInterpOfOne);
		_impl->Setup(*this, &a, a_TT.Get(), a_TF.Get(), nullptr, nullptr, cellInterpOfOne, faceInterpOfOne);
	}

	template <typename SparseMatrixType>
	void Solver::SetupSymmetric(const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FF,
	                            const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
	{
		CallScope scope(*this);
		_impl->CheckParameters(*this, true);
		InputMatrix a_TT(A_TT, "A_TT"), a_TF(A_TF, "A_TF"), a_FF(A_FF, "A_FF");
		_impl->CheckInputs(*this, nullptr, a_TT.Get(), a_TF.Get(), nullptr, &a_FF.Get(), cellInterpOfOne, faceInterpOfOne);
		_impl->Setup(*this, nullptr, a_TT.Get(), a_TF.Get(), nullptr, &a_FF.Get(), cellInterpOfOne, faceInterpOfOne);
	}

	template <typename SparseMatrixType>
	void Solver::SetupNonsymmetricCondensed(const SparseMatrixType& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FT, const SparseMatrixType& A_FF,
	                                        const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
	{
		CallScope scope(*this);
		_impl->CheckParameters(*this, false);
		InputMatrix a(A, "A"), a_TT(A_TT, "A_TT"), a_TF(A_TF, "A_TF"), a_FT(A_FT, "A_FT"), a_FF(A_FF, "A_FF");
		_impl->CheckInputs(*this, &a.Get(), a_TT.Get(), a_TF.Get(), &a_FT.Get(), &a_FF.Get(), cellInterpOfOne, faceInterpOfOne);
		_impl->Setup(*this, &a, a_TT.Get(), a_TF.Get(), &a_FT.Get(), &a_FF.Get(), cellInterpOfOne, faceInterpOfOne);
	}

	template <typename SparseMatrixType>
	void Solver::SetupNonsymmetricCondensed(const SparseMatrixType& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FT,
	                                        const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
	{
		CallScope scope(*this);
		_impl->CheckParameters(*this, false);
		InputMatrix a(A, "A"), a_TT(A_TT, "A_TT"), a_TF(A_TF, "A_TF"), a_FT(A_FT, "A_FT");
		_impl->CheckInputs(*this, &a.Get(), a_TT.Get(), a_TF.Get(), &a_FT.Get(), nullptr, cellInterpOfOne, faceInterpOfOne);
		_impl->Setup(*this, &a, a_TT.Get(), a_TF.Get(), &a_FT.Get(), nullptr, cellInterpOfOne, faceInterpOfOne);
	}

	template <typename SparseMatrixType>
	void Solver::SetupNonsymmetric(const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FT, const SparseMatrixType& A_FF,
	                               const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
	{
		CallScope scope(*this);
		_impl->CheckParameters(*this, false);
		InputMatrix a_TT(A_TT, "A_TT"), a_TF(A_TF, "A_TF"), a_FT(A_FT, "A_FT"), a_FF(A_FF, "A_FF");
		_impl->CheckInputs(*this, nullptr, a_TT.Get(), a_TF.Get(), &a_FT.Get(), &a_FF.Get(), cellInterpOfOne, faceInterpOfOne);
		_impl->Setup(*this, nullptr, a_TT.Get(), a_TF.Get(), &a_FT.Get(), &a_FF.Get(), cellInterpOfOne, faceInterpOfOne);
	}

	Result Solver::Solve(const Eigen::Ref<const Eigen::VectorXd>& b, Eigen::Ref<Eigen::VectorXd> x)
	{
		_impl->CheckSetUp("Solve");
		_impl->CheckVectorSize(b.rows(), "b");
		_impl->CheckVectorSize(x.rows(), "x");
		CallScope scope(*this);
		auto start = chrono::steady_clock::now();

		Result result;
		Vector rhs = b;
		if (rhs.norm() == 0)
		{
			x.setZero();
			result.Converged = true;
			return result;
		}

		Vector solution = x;
		bool xEquals0 = (solution.array() == 0).all(); // skips the products by the initial guess, as the program does
		IterativeSolver* solver = _impl->ActiveSolver();
		solver->Tolerance = Tolerance;
		solver->MaxIterations = MaxIterations;
		solver->Solve(rhs, solution, xEquals0);
		x = solution;

		result.Iterations = solver->IterationCount;
		result.RelativeResidual = solver->LastIterationResult.NormalizedResidualNorm;
		result.Converged = result.RelativeResidual < Tolerance;
		if (Verbosity == 1)
			cout << "fhhos4: " << result.Iterations << " iterations, relative residual " << result.RelativeResidual << (result.Converged ? "" : " (not converged)") << ", " << Seconds(start) << " s" << endl;
		return result;
	}

	Result Solver::Solve(const Eigen::Ref<const Eigen::VectorXd>& b_T, const Eigen::Ref<const Eigen::VectorXd>& b_F,
	                     Eigen::Ref<Eigen::VectorXd> x_T, Eigen::Ref<Eigen::VectorXd> x_F)
	{
		_impl->CheckSetUp("Solve");
		if (!_impl->Condensation)
			Utils::FatalError("fhhos4::Solver::Solve(b_T, b_F, x_T, x_F) solves the uncondensed system: it needs a setup from its blocks (SetupSymmetric() or SetupNonsymmetric()). After SetupSymmetricCondensed() or SetupNonsymmetricCondensed(), solve the condensed system with Solve(b, x).");
		_impl->CheckCellVectorSize(b_T.rows(), "b_T");
		_impl->CheckCellVectorSize(x_T.rows(), "x_T");
		_impl->CheckVectorSize(b_F.rows(), "b_F");
		_impl->CheckVectorSize(x_F.rows(), "x_F");
		CallScope scope(*this);

		const StaticCondensation& condensation = *_impl->Condensation;
		Vector invA_TT_b_T = condensation.SolveCells(b_T);
		Vector b = condensation.CondensedRHS(b_F, invA_TT_b_T);
		Result result = Solve(b, x_F);
		x_T = condensation.CellUnknowns(invA_TT_b_T, x_F);
		return result;
	}

	const SparseMatrix& Solver::CondensedMatrix() const
	{
		_impl->CheckSetUp("CondensedMatrix");
		return *_impl->A;
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

	// The supported matrices
#define FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex) const Eigen::SparseMatrix<double, StorageOrder, StorageIndex>&
#define FHHOS4_AMG_INSTANTIATE_SETUP(StorageOrder, StorageIndex) \
	template void Solver::SetupSymmetricCondensed<Eigen::SparseMatrix<double, StorageOrder, StorageIndex>>( \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&); \
	template void Solver::SetupSymmetricCondensed<Eigen::SparseMatrix<double, StorageOrder, StorageIndex>>( \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&); \
	template void Solver::SetupSymmetric<Eigen::SparseMatrix<double, StorageOrder, StorageIndex>>( \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&); \
	template void Solver::SetupNonsymmetricCondensed<Eigen::SparseMatrix<double, StorageOrder, StorageIndex>>( \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&); \
	template void Solver::SetupNonsymmetricCondensed<Eigen::SparseMatrix<double, StorageOrder, StorageIndex>>( \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&); \
	template void Solver::SetupNonsymmetric<Eigen::SparseMatrix<double, StorageOrder, StorageIndex>>( \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&);

	FHHOS4_AMG_INSTANTIATE_SETUP(Eigen::ColMajor, int)
	FHHOS4_AMG_INSTANTIATE_SETUP(Eigen::RowMajor, int)
	FHHOS4_AMG_INSTANTIATE_SETUP(Eigen::ColMajor, long)
	FHHOS4_AMG_INSTANTIATE_SETUP(Eigen::RowMajor, long)
}
