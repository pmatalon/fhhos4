#include "fhhos4"
#include "LibraryUAMG.h"

class LibraryUAMG : public IterativeSolver
{
private:
	fhhos4::Solver _solver;
public:
	// The parameters of the library that correspond to the arguments of the program
	LibraryUAMG(const ProgramArguments& args, int dim, int faceDegree, int cellDegree, bool useFCG)
	{
		_solver.Dimension = dim;
		_solver.FaceDegree = faceDegree;
		_solver.CellDegree = cellDegree;
		_solver.Krylov = useFCG ? "fcg" : "none";
		_solver.Tolerance = args.Solver.Tolerance;
		_solver.MaxIterations = args.Solver.MaxIterations;
		const MultigridArguments& mg = args.Solver.MG;
		_solver.Cycle = mg.CycleLetter;
		_solver.PreSmoothingIterations = mg.PreSmoothingIterations;
		_solver.PostSmoothingIterations = mg.PostSmoothingIterations;
		_solver.PreSmoother = mg.PreSmootherCode;
		_solver.PostSmoother = mg.PostSmootherCode;
		_solver.CoarseningFactor = mg.CoarseningFactor;
		_solver.CoarseMatrixMaxSize = mg.MatrixMaxSizeForCoarsestLevel;
		_solver.CoarseSolver = mg.CoarseSolverCode;
		_solver.Threads = 0; // the setting of the program (-threads)
		_solver.Verbosity = args.Solver.PrintIterationResults ? 2 : 1;
	}

	void Serialize(ostream& os) const override
	{
		const fhhos4::Solver& s = _solver;
		os << "U-AMG through the library fhhos4_AMG, Krylov method: " << s.Krylov << endl;
		os << "\t" << "Cycle                   : " << s.Cycle << "(" << s.PreSmoothingIterations << "," << s.PostSmoothingIterations << ")" << endl;
		os << "\t" << "Smoothers               : " << s.PreSmoother << ", " << s.PostSmoother << endl;
		os << "\t" << "Coarsening factor       : " << s.CoarseningFactor << endl;
		os << "\t" << "Coarse solver           : " << s.CoarseSolver << " (matrix size <= " << s.CoarseMatrixMaxSize << ")";
	}

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
		IterativeSolver::Setup(A);
		// Without A_F_F, which the library computes from A: the setup of the codes that have A (e.g. HArDCore3D)
		_solver.Setup(A, A_T_T, A_T_F, cellInterpOfOne, faceInterpOfOne);
	}

	void Solve(const Vector& b, Vector& x, bool xEquals0, bool computeResidual, bool computeAx) override
	{
		if (computeResidual || computeAx)
			Utils::FatalError("LibraryUAMG: computeResidual and computeAx are not managed.");

		fhhos4::Result result = _solver.Solve(b, x);
		this->IterationCount = result.Iterations;
		this->LastIterationResult.NormalizedResidualNorm = result.RelativeResidual;
	}
};

IterativeSolver* CreateLibraryUAMG(const ProgramArguments& args, int dim, int faceDegree, int cellDegree, bool useFCG)
{
	return new LibraryUAMG(args, dim, faceDegree, cellDegree, useFCG);
}
