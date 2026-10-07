// The shared library fhhos4_AMG: the algebraic multigrid of fhhos4 for hybrid discretizations (U-AMG), for the codes
// that discretize the problem themselves (include/fhhos4). Compiled with fhhos4's flags and -fvisibility=hidden: only
// the classes of the public header are exported. The multigrid is built as in the program (SolverFactory,
// ApplyProgramArgumentDefaults), from the same code (AlgebraicSolverFactory, ApplyUncondensedAMGDefaults).
#include <chrono>
#include <limits>
#include "fhhos4"
#include "ProgramArgumentsDefaults.h"
#include "Solver/AlgebraicSolverFactory.h"
#include "Solver/Krylov/PreconditionedBiCGSTAB.h"
#include "Utils/SparseMatrixOps.h"

namespace fhhos4
{
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
				if (solver.Threads > 0)
				{
					_savedThreads = omp_get_max_threads();
					omp_set_num_threads(solver.Threads);
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

		// Copy of a matrix of the calling code in fhhos4's format (row-major, fhhos4's index type, compressed)
		template <typename SparseMatrixType>
		SparseMatrix Copy(const SparseMatrixType& M, const string& name)
		{
			if (M.nonZeros() > (Eigen::Index)numeric_limits<SparseMatrixIndex>::max())
				Utils::FatalError("fhhos4::Solver: " + name + " has too many non-zeros for the index type of fhhos4 (compile fhhos4 with -DSMALL_INDEX=OFF).");
			SparseMatrix copy = M;
			copy.makeCompressed();
			return copy;
		}

		// A_FT A_TT^{-1} A_TF, from the lower part of A_TT
		SparseMatrix CellElimination(const SparseMatrix& A_TT, const SparseMatrix& A_TF, int cellBlockSize)
		{
			SparseMatrix invA_TT = Utils::InvertBlockDiagMatrix(SparseMatrixOps::FullFromLower(A_TT), cellBlockSize);
			return SparseMatrixOps::Multiply(SparseMatrixOps::Transpose(A_TF), SparseMatrixOps::Multiply(invA_TT, A_TF));
		}

		// A = A_FF - A_FT A_TT^{-1} A_TF, from the lower parts of A_TT and A_FF
		SparseMatrix CondensedMatrix(const SparseMatrix& A_TT, const SparseMatrix& A_TF, const SparseMatrix& A_FF, int cellBlockSize)
		{
			SparseMatrix A = SparseMatrixOps::FullFromLower(A_FF) - CellElimination(A_TT, A_TF, cellBlockSize);
			A.makeCompressed();
			return A;
		}

		// A_FF = A + A_FT A_TT^{-1} A_TF, from the lower part of A_TT
		SparseMatrix FaceBlock(const SparseMatrix& A, const SparseMatrix& A_TT, const SparseMatrix& A_TF, int cellBlockSize)
		{
			SparseMatrix A_FF = A + CellElimination(A_TT, A_TF, cellBlockSize);
			A_FF.makeCompressed();
			return A_FF;
		}
	}

	struct Solver::Impl
	{
		// Set by Setup()
		int CellDegree = 0;
		int CellBlockSize = 0;
		int FaceBlockSize = 0;
		SparseMatrix A; // the solvers keep a pointer to it
		unique_ptr<UncondensedAMG> Multigrid;
		unique_ptr<IterativeSolver> Krylov; // preconditioned by Multigrid; null if Solver::Krylov is "none"
		bool IsSetUp = false;

		IterativeSolver* ActiveSolver() const
		{
			if (Krylov)
				return Krylov.get();
			return Multigrid.get();
		}

		// Checks the parameters, and computes the block sizes
		void CheckParameters(const Solver& s)
		{
			if (s.Dimension != 2 && s.Dimension != 3)
				Utils::FatalError("fhhos4::Solver: set Dimension (2 or 3) before Setup(), got " + to_string(s.Dimension) + ".");
			if (s.FaceDegree < 0)
				Utils::FatalError("fhhos4::Solver: set FaceDegree (>= 0) before Setup(), got " + to_string(s.FaceDegree) + ".");
			if (s.CellDegree < -1)
				Utils::FatalError("fhhos4::Solver: CellDegree must be >= 0 (or -1: FaceDegree), got " + to_string(s.CellDegree) + ".");
			if (s.Cycle != 'V' && s.Cycle != 'W' && s.Cycle != 'K')
				Utils::FatalError(string("fhhos4::Solver: unknown Cycle '") + s.Cycle + "' (V, W or K).");
			if (s.Krylov != "fcg" && s.Krylov != "bicgstab" && s.Krylov != "none")
				Utils::FatalError("fhhos4::Solver: unknown Krylov method '" + s.Krylov + "' (fcg, bicgstab or none).");
			if (s.Krylov == "bicgstab" && s.Cycle == 'K')
				Utils::FatalError("fhhos4::Solver: BiCGSTAB is not flexible: it needs a linear preconditioner, a V- or W-cycle, not the K-cycle (non-linear: it solves the coarse levels with FCG iterations). Set Cycle to 'V' or 'W', or Krylov to \"fcg\".");
			CellDegree = s.CellDegree == -1 ? s.FaceDegree : s.CellDegree;
			CellBlockSize = Utils::Binomial(CellDegree + s.Dimension, CellDegree);
			FaceBlockSize = Utils::Binomial(s.FaceDegree + s.Dimension - 1, s.FaceDegree);
		}

		// The sizes and structure of the inputs, against the dimension and the degrees. A, A_FF: null if not given.
		void CheckInputs(const Solver& s, const SparseMatrix* A, const SparseMatrix& A_TT, const SparseMatrix& A_TF, const SparseMatrix* A_FF, const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne) const
		{
			string cellDoFs = to_string(CellBlockSize) + " DoFs per cell (CellDegree " + to_string(CellDegree) + ", Dimension " + to_string(s.Dimension) + ")";
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
			if (A_FF && (A_FF->rows() != nFaceDoFs || A_FF->cols() != nFaceDoFs))
				Utils::FatalError("fhhos4::Solver: A_FF is " + size(*A_FF) + ", expected " + to_string(nFaceDoFs) + " x " + to_string(nFaceDoFs) + " (the columns of A_TF).");
			if (A && (A->rows() != nFaceDoFs || A->cols() != nFaceDoFs))
				Utils::FatalError("fhhos4::Solver: A is " + size(*A) + ", expected " + to_string(nFaceDoFs) + " x " + to_string(nFaceDoFs) + " (the columns of A_TF).");
			if (cellInterpOfOne.rows() != nCellDoFs || faceInterpOfOne.rows() != nFaceDoFs)
				Utils::FatalError("fhhos4::Solver: the interpolations of 1 on the cells and faces must have the sizes of the rows (" + to_string(nCellDoFs) + ") and columns (" + to_string(nFaceDoFs) + ") of A_TF, got " + to_string(cellInterpOfOne.rows()) + " and " + to_string(faceInterpOfOne.rows()) + ".");

			// A_TT is block diagonal
			for (BigNumber i = 0; i < A_TT.rows(); i++)
			{
				for (SparseMatrix::InnerIterator it(A_TT, i); it; ++it)
				{
					if (it.col() / CellBlockSize != i / CellBlockSize && it.value() != 0)
						Utils::FatalError("fhhos4::Solver: A_TT has a non-zero (" + to_string(i) + ", " + to_string(it.col()) + ") outside the diagonal blocks of the cells, with " + cellDoFs + ".");
				}
			}
			// The interpolations of 1: non-zero on the first DoF of each cell or face, zero on the others
			CheckInterpOfOne(cellInterpOfOne, CellBlockSize, "cell", cellDoFs);
			CheckInterpOfOne(faceInterpOfOne, FaceBlockSize, "face", faceDoFs);
		}

		static void CheckInterpOfOne(const Eigen::Ref<const Eigen::VectorXd>& interpOfOne, int blockSize, const string& entity, const string& blockDoFs)
		{
			for (Eigen::Index block = 0; block * blockSize < interpOfOne.rows(); block++)
			{
				double c = interpOfOne[block * blockSize];
				if (!std::isfinite(c) || c == 0)
					Utils::FatalError("fhhos4::Solver: the interpolation of 1 on the " + entity + " bases is not finite and non-zero on the first DoF of the " + entity + " " + to_string(block) + ", with " + blockDoFs + ".");
				for (int i = 1; i < blockSize; i++)
				{
					if (!(abs(interpOfOne[block * blockSize + i]) <= 1e-10 * abs(c)))
						Utils::FatalError("fhhos4::Solver: the interpolation of 1 on the " + entity + " bases is not zero beyond the first DoF of the " + entity + " " + to_string(block) + ", with " + blockDoFs + ". The bases must be hierarchical, with a constant first function.");
				}
			}
		}

		// The arguments of the program fhhos4 that give the same multigrid
		ProgramArguments Arguments(const Solver& s) const
		{
			ProgramArguments args;
			args.Solver.Tolerance = s.Tolerance;
			args.Solver.MaxIterations = s.MaxIterations;
			args.Solver.PrintIterationResults = s.Verbosity >= 2;
			MultigridArguments& mg = args.Solver.MG;
			mg.CycleLetter = s.Cycle;
			mg.WLoops = s.Cycle == 'W' ? 2 : 1;
			mg.PreSmoothingIterations = s.PreSmoothingIterations;
			mg.PostSmoothingIterations = s.PostSmoothingIterations;
			mg.PreSmootherCode = s.PreSmoother;
			mg.PostSmootherCode = s.PostSmoother;
			mg.CoarseningFactor = s.CoarseningFactor;
			mg.MatrixMaxSizeForCoarsestLevel = s.CoarseMatrixMaxSize;
			mg.CoarseSolverCode = s.CoarseSolver;
			bool defaultCoarseOperator = true, defaultCycle = false, defaultHPCoarseningStgy = true;
			ApplyUncondensedAMGDefaults(mg, s.FaceDegree, defaultCoarseOperator, defaultCycle, defaultHPCoarseningStgy);
			return args;
		}

		// A, A_FF: moved; null to compute them from the other matrices (not both)
		void Setup(const Solver& s, SparseMatrix* A, SparseMatrix&& A_TT, SparseMatrix&& A_TF, SparseMatrix* A_FF, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne)
		{
			auto start = chrono::steady_clock::now();
			IsSetUp = false;
			Krylov.reset();
			Multigrid.reset();

			ProgramArguments args = Arguments(s);
			ExportModule out;
			Multigrid.reset(AlgebraicSolverFactory::CreateUncondensedAMG(args, s.Dimension, s.FaceDegree, CellBlockSize, FaceBlockSize, FaceBlockSize, out, AlgebraicSolverFactory::AlgebraicCoarseSolvers(s.Dimension, out)));

			SparseMatrix faceBlock; // empty if neither given nor needed (the algorithm of the paper only uses A_TT and A_TF)
			if (A_FF)
				faceBlock = std::move(*A_FF);
			if (A)
				this->A = std::move(*A);
			else
				this->A = CondensedMatrix(A_TT, A_TF, faceBlock, CellBlockSize);
			if (!A_FF && Multigrid->A_F_FNeeded())
				faceBlock = FaceBlock(this->A, A_TT, A_TF, CellBlockSize);

			// The Krylov method preconditioned by one cycle (SolverPreconditioner), as SolverFactory (-s fcguamg)
			if (s.Krylov == "fcg")
			{
				auto fcg = make_unique<FlexibleConjugateGradient>(1);
				fcg->Precond = SolverPreconditioner(Multigrid.get());
				Krylov = std::move(fcg);
			}
			else if (s.Krylov == "bicgstab")
			{
				auto bicgstab = make_unique<PreconditionedBiCGSTAB>();
				bicgstab->Precond = SolverPreconditioner(Multigrid.get());
				Krylov = std::move(bicgstab);
			}
			ActiveSolver()->PrintIterationResults = args.Solver.PrintIterationResults;

			// The blocks are only read during the setup (UncondensedAMG::Setup())
			ActiveSolver()->Setup(this->A, A_TT, A_TF, faceBlock, cellInterpOfOne, faceInterpOfOne);
			IsSetUp = true;

			if (s.Verbosity == 1)
				cout << "fhhos4_AMG: setup of " << Multigrid->NumberOfLevels() << " levels in " << Seconds(start) << " s" << endl;
		}

		void CheckSetUp(const char* function) const
		{
			if (!IsSetUp)
				Utils::FatalError(string("fhhos4::Solver::") + function + "(): Setup() must be called first (and succeed).");
		}

		void CheckVectorSize(Eigen::Index size, const char* vector) const
		{
			if (size != A.rows())
				Utils::FatalError(string("fhhos4::Solver: the vector ") + vector + " has " + to_string(size) + " rows, expected " + to_string(A.rows()) + " (the faces DoFs).");
		}
	};

	Solver::Solver(const EigenConfiguration& host)
	{
		EigenConfiguration library = HostEigenConfiguration(); // the configuration of the library's compilation
		auto toString = [](const EigenConfiguration& c)
		{
			return to_string(c.World) + "." + to_string(c.Major) + "." + to_string(c.Minor) + ", max alignment " + to_string(c.MaxAlignBytes) + " bytes, index of " + to_string(c.IndexSize) + " bytes";
		};
		if (host.World != library.World || host.Major != library.Major || host.Minor != library.Minor || host.MaxAlignBytes != library.MaxAlignBytes || host.IndexSize != library.IndexSize)
			throw Error("fhhos4::Solver: the calling code and the library fhhos4_AMG are compiled with incompatible configurations of Eigen (calling code: Eigen " + toString(host) + "; library: Eigen " + toString(library) + "). Use the same version of Eigen, and the same -march (it sets the alignment).");
		_impl = make_unique<Impl>();
	}

	Solver::~Solver() = default;
	Solver::Solver(Solver&&) noexcept = default;
	Solver& Solver::operator=(Solver&&) noexcept = default;

	template <typename SparseMatrixType>
	void Solver::Setup(const SparseMatrixType& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FF,
	                   const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
	{
		CallScope scope(*this);
		_impl->CheckParameters(*this);
		SparseMatrix copyA = Copy(A, "A");
		SparseMatrix copyA_TT = Copy(A_TT, "A_TT");
		SparseMatrix copyA_TF = Copy(A_TF, "A_TF");
		SparseMatrix copyA_FF = Copy(A_FF, "A_FF");
		_impl->CheckInputs(*this, &copyA, copyA_TT, copyA_TF, &copyA_FF, cellInterpOfOne, faceInterpOfOne);
		_impl->Setup(*this, &copyA, std::move(copyA_TT), std::move(copyA_TF), &copyA_FF, cellInterpOfOne, faceInterpOfOne);
	}

	template <typename SparseMatrixType>
	void Solver::Setup(const SparseMatrixType& A, const SparseMatrixType& A_TT, const SparseMatrixType& A_TF,
	                   const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
	{
		CallScope scope(*this);
		_impl->CheckParameters(*this);
		SparseMatrix copyA = Copy(A, "A");
		SparseMatrix copyA_TT = Copy(A_TT, "A_TT");
		SparseMatrix copyA_TF = Copy(A_TF, "A_TF");
		_impl->CheckInputs(*this, &copyA, copyA_TT, copyA_TF, nullptr, cellInterpOfOne, faceInterpOfOne);
		_impl->Setup(*this, &copyA, std::move(copyA_TT), std::move(copyA_TF), nullptr, cellInterpOfOne, faceInterpOfOne);
	}

	template <typename SparseMatrixType>
	void Solver::SetupFromBlocks(const SparseMatrixType& A_TT, const SparseMatrixType& A_TF, const SparseMatrixType& A_FF,
	                             const Eigen::Ref<const Eigen::VectorXd>& cellInterpOfOne, const Eigen::Ref<const Eigen::VectorXd>& faceInterpOfOne)
	{
		CallScope scope(*this);
		_impl->CheckParameters(*this);
		SparseMatrix copyA_TT = Copy(A_TT, "A_TT");
		SparseMatrix copyA_TF = Copy(A_TF, "A_TF");
		SparseMatrix copyA_FF = Copy(A_FF, "A_FF");
		_impl->CheckInputs(*this, nullptr, copyA_TT, copyA_TF, &copyA_FF, cellInterpOfOne, faceInterpOfOne);
		_impl->Setup(*this, nullptr, std::move(copyA_TT), std::move(copyA_TF), &copyA_FF, cellInterpOfOne, faceInterpOfOne);
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
			cout << "fhhos4_AMG: " << result.Iterations << " iterations, relative residual " << result.RelativeResidual << (result.Converged ? "" : " (not converged)") << ", " << Seconds(start) << " s" << endl;
		return result;
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
	template void Solver::Setup<Eigen::SparseMatrix<double, StorageOrder, StorageIndex>>( \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&); \
	template void Solver::Setup<Eigen::SparseMatrix<double, StorageOrder, StorageIndex>>( \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&); \
	template void Solver::SetupFromBlocks<Eigen::SparseMatrix<double, StorageOrder, StorageIndex>>( \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		FHHOS4_AMG_MATRIX(StorageOrder, StorageIndex), \
		const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&);

	FHHOS4_AMG_INSTANTIATE_SETUP(Eigen::ColMajor, int)
	FHHOS4_AMG_INSTANTIATE_SETUP(Eigen::RowMajor, int)
	FHHOS4_AMG_INSTANTIATE_SETUP(Eigen::ColMajor, long)
	FHHOS4_AMG_INSTANTIATE_SETUP(Eigen::RowMajor, long)
}
