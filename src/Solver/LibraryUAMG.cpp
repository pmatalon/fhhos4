#include "fhhos4.h"
#include "LibraryUAMG.h"
#include "LibraryParameters.h"

class LibraryUAMGImpl : public LibraryUAMG
{
private:
	fhhos4::Solver _solver;
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
		if (!krylov)
			_solver.Krylov = fhhos4::Krylov::None;
		else
			_solver.Krylov = *krylov == KrylovMethod::GCR ? fhhos4::Krylov::GCR : fhhos4::Krylov::FCG;
		_solver.Tolerance = args.Solver.Tolerance;
		_solver.MaxIterations = args.Solver.MaxIterations;
		const MultigridArguments& mg = args.Solver.MG;
		_solver.Cycle = LibraryParameters::LibraryCycle(mg.Cycle);
		_solver.PCoarsening = LibraryParameters::LibraryPCoarsening(mg.P_CS);
		_solver.PreSmoothingIterations = mg.PreSmoothingIterations;
		_solver.PostSmoothingIterations = mg.PostSmoothingIterations;
		_solver.PreSmoother = LibraryParameters::LibrarySmoother(mg.PreSmoother);
		_solver.PostSmoother = LibraryParameters::LibrarySmoother(mg.PostSmoother);
		_solver.CoarseningFactor = mg.CoarseningFactor;
		_solver.CoarseMatrixMaxSize = mg.MatrixMaxSizeForCoarsestLevel;
		_solver.CoarseSolver = LibraryParameters::LibraryCoarseSolver(mg.CoarseSolverCode);
		_solver.Threads = Parallelism::NThreads(); // the setting of the program for its solvers (-threads, SolverThreads)
		_solver.Verbosity = args.Solver.PrintIterationResults ? 2 : 1;
	}

	void Serialize(ostream& os) const override
	{
		os << "U-AMG through the library fhhos4, Krylov method: " << (!_krylov ? "none" : *_krylov == KrylovMethod::GCR ? "gcr" : "fcg") << endl;
		os << "\t" << "Cycle                   : " << (_mg.Cycle == CycleType::K ? "K" : _mg.Cycle == CycleType::W ? "W" : "V") << "(" << _mg.PreSmoothingIterations << "," << _mg.PostSmoothingIterations << ")" << endl;
		os << "\t" << "Smoothers               : " << SmootherFactory::Code(_mg.PreSmoother) << ", " << SmootherFactory::Code(_mg.PostSmoother) << endl;
		os << "\t" << "Coarsening factor       : " << _mg.CoarseningFactor << endl;
		os << "\t" << "Coarse solver           : " << _mg.CoarseSolverCode << " (matrix size <= " << _mg.MatrixMaxSizeForCoarsestLevel << ")";
	}

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
		IterativeSolver::Setup(A);
		// Without A_F_F, which the library computes from A: the setup of the codes that have A (e.g. HArDCore3D)
		if (A_F_T.rows() > 0)
			_solver.SetupNonsymmetricCondensed(A, A_T_T, A_T_F, A_F_T, cellInterpOfOne, faceInterpOfOne);
		else
			_solver.SetupSymmetricCondensed(A, A_T_T, A_T_F, cellInterpOfOne, faceInterpOfOne);
	}

	void SetupFromBlocks(const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
		if (A_F_T.rows() > 0)
			_solver.SetupNonsymmetric(A_T_T, A_T_F, A_F_T, A_F_F, cellInterpOfOne, faceInterpOfOne);
		else
			_solver.SetupSymmetric(A_T_T, A_T_F, A_F_F, cellInterpOfOne, faceInterpOfOne);
		IterativeSolver::Setup(_solver.CondensedMatrix());
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
