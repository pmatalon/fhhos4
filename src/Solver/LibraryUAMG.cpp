#include "fhhos4.h"
#include "LibraryUAMG.h"
#include "LibraryParameters.h"
#include "../Utils/BorrowedSparseMatrix.h"
#include "../Utils/Communicator.h"

class LibraryUAMGImpl : public LibraryUAMG
{
private:
	fhhos4::Solver _solver;
	BorrowedSparseMatrix _condensedMatrix; // the library's: the matrix of the solver
	MultigridArguments _mg;         // printed by Serialize(), with the codes of the program
	optional<KrylovMethod> _krylov;
public:
	// The parameters of the library that correspond to the arguments of the program
	LibraryUAMGImpl(const ProgramArguments& args, int dim, int faceDegree, int cellDegree, optional<KrylovMethod> krylov)
		: _mg(args.Solver.MG), _krylov(krylov)
	{
		_solver.Dimension = dim;
		_solver.FaceDegree = faceDegree;
		_solver.CellDegree = cellDegree;
		_solver.PolynomialSpace = LibraryParameters::LibraryPolynomialSpace(args.Discretization.UsePolynomialSpaceQ);
		if (!krylov)
			_solver.Krylov = fhhos4::Krylov::None;
		else
			_solver.Krylov = *krylov == KrylovMethod::GCR ? fhhos4::Krylov::GCR : fhhos4::Krylov::FCG;
		_solver.Tolerance = args.Solver.Tolerance;
		_solver.MaxIterations = args.Solver.MaxIterations;
		const MultigridArguments& mg = args.Solver.MG;
		_solver.Cycle = LibraryParameters::LibraryCycle(mg.Cycle);
		_solver.PCoarsening = LibraryParameters::LibraryPCoarsening(mg.P_CS);
		_solver.Prolongation = LibraryParameters::LibraryProlongation(mg.UAMGCoarseningProlong);
		_solver.PreSmoothingIterations = mg.PreSmoothingIterations;
		_solver.PostSmoothingIterations = mg.PostSmoothingIterations;
		_solver.PreSmoother = LibraryParameters::LibrarySmoother(mg.PreSmoother);
		_solver.PostSmoother = LibraryParameters::LibrarySmoother(mg.PostSmoother);
		_solver.CoarseningFactor = mg.CoarseningFactor;
		_solver.CoarseMatrixMaxSize = mg.MatrixMaxSizeForCoarsestLevel;
		_solver.CoarseSolver = LibraryParameters::LibraryCoarseSolver(mg.CoarseSolverCode);
		_solver.MixedPrecision = mg.UAMGMixedPrecision; // the program's default: double (the papers)
		_solver.BlockCSRStorage = mg.UAMGBlockCSRStorage;
		_solver.Threads = Parallelism::NThreads(); // the setting of the program for its solvers (-threads, SolverThreads)
		_solver.Verbosity = args.Solver.PrintIterationResults ? 2 : 1;
#if FHHOS4_MPI
		_solver.Comm = MPI_COMM_SELF; // the system of the program on one process; SetupFromDistributedBlocks(): all of them
#endif
	}

	void Serialize(ostream& os) const override
	{
		os << "U-AMG through the library fhhos4, Krylov method: " << (!_krylov ? "none" : *_krylov == KrylovMethod::GCR ? "gcr" : "fcg") << endl;
		os << "\t" << "Cycle                   : " << (_mg.Cycle == CycleType::K ? "K" : _mg.Cycle == CycleType::W ? "W" : "V") << "(" << _mg.PreSmoothingIterations << "," << _mg.PostSmoothingIterations << ")" << endl;
		os << "\t" << "Smoothers               : " << SmootherFactory::Code(_mg.PreSmoother) << ", " << SmootherFactory::Code(_mg.PostSmoother) << endl;
		os << "\t" << "Coarsening prolongation : " << (unsigned)_mg.UAMGCoarseningProlong << endl;
		os << "\t" << "Coarsening factor       : " << _mg.CoarseningFactor << endl;
		os << "\t" << "Coarse solver           : " << _mg.CoarseSolverCode << " (matrix size <= " << _mg.MatrixMaxSizeForCoarsestLevel << ")";
	}

	// The program solves the condensed system (Solve(b)): the setup from the blocks, A computed by the library (equal to
	// the program's up to rounding)
	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
		SetupFromBlocks(A_T_T, A_T_F, A_F_T, A_F_F, cellInterpOfOne, faceInterpOfOne);
	}

	void SetupFromBlocks(const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
		if (A_F_T.rows() > 0)
			_solver.SetupNonsymmetric(A_T_T, A_T_F, A_F_T, A_F_F, cellInterpOfOne, faceInterpOfOne);
		else
			_solver.SetupSymmetric(A_T_T, A_T_F, A_F_F, cellInterpOfOne, faceInterpOfOne);
		fhhos4::SparseMatrixMap A = _solver.CondensedMatrix();
		_condensedMatrix.Borrow(A.rows(), A.cols(), A.outerIndexPtr(), A.innerIndexPtr(), A.valuePtr());
		IterativeSolver::Setup(_condensedMatrix);
	}

	void SetupFromDistributedBlocks(const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F, const vector<int64_t>& faceGhosts, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
#ifdef MPI_ENABLED
		_solver.Comm = Communicator::Current().Comm();
#endif
		_solver.SetupSymmetric(A_T_T, A_T_F, A_F_F, cellInterpOfOne, faceInterpOfOne, fhhos4::GhostNumbers(faceGhosts.size(), faceGhosts.data()));
	}

	void Solve(const Vector& b, Vector& x, bool xEquals0, bool computeResidual, bool computeAx) override
	{
		if (computeResidual || computeAx)
			Utils::FatalError("LibraryUAMG: computeResidual and computeAx are not managed.");

		fhhos4::Result result = _solver.Solve(b, x);
		this->IterationCount = result.Iterations;
		this->LastIterationResult.NormalizedResidualNorm = result.RelativeResidual;
	}

	void Solve(const Vector& b_T, const Vector& b_F, Vector& x_T, Vector& x_F) override
	{
		x_T.resize(b_T.rows());
		fhhos4::Result result = _solver.Solve(b_T, b_F, x_T, x_F);
		this->IterationCount = result.Iterations;
		this->LastIterationResult.NormalizedResidualNorm = result.RelativeResidual;
	}
};

LibraryUAMG* CreateLibraryUAMG(const ProgramArguments& args, int dim, int faceDegree, int cellDegree, optional<KrylovMethod> krylov)
{
	return new LibraryUAMGImpl(args, dim, faceDegree, cellDegree, krylov);
}
