#pragma once
#include <Eigen/IterativeLinearSolvers>
#include "../IterativeSolver.h"
#include "../Preconditioner.h"
using namespace std;

// Eigen's BiCGSTAB preconditioned by a solver of fhhos4 (e.g. one cycle of a multigrid). BiCGSTAB is not flexible: the
// preconditioner must be a fixed linear operator (e.g. a V- or W-cycle, not a K-cycle). EigenBiCGSTAB has the diagonal
// preconditioner.
class PreconditionedBiCGSTAB : public IterativeSolver
{
private:
	// The preconditioner in the form that Eigen's iterative solvers expect
	class EigenPreconditioner
	{
	public:
		SolverPreconditioner* Precond = nullptr;

		template <typename MatrixType>
		EigenPreconditioner& analyzePattern(const MatrixType&) { return *this; }
		template <typename MatrixType>
		EigenPreconditioner& factorize(const MatrixType&) { return *this; }
		template <typename MatrixType>
		EigenPreconditioner& compute(const MatrixType&) { return *this; }

		template <typename Rhs>
		Vector solve(const Rhs& r) const
		{
			return Precond->Apply(Vector(r));
		}

		Eigen::ComputationInfo info() const
		{
			return Eigen::Success;
		}
	};

	Eigen::BiCGSTAB<SparseMatrix, EigenPreconditioner> _solver;

public:
	SolverPreconditioner Precond;

	void Serialize(ostream& os) const override
	{
		os << "BiCGSTAB (Eigen library), preconditioner: " << Precond;
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		this->Precond.Setup(A);
		EndSetup(A);
	}

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
		IterativeSolver::Setup(A);
		this->Precond.Setup(A, A_T_T, A_T_F, A_F_F, cellInterpOfOne, faceInterpOfOne);
		EndSetup(A);
	}

	void Solve(const Vector& b, Vector& x, bool xEquals0, bool computeResidual, bool computeAx) override
	{
		if (computeResidual || computeAx)
			Utils::FatalError("PreconditionedBiCGSTAB: computeResidual and computeAx are not managed.");
		_solver.setMaxIterations(this->MaxIterations);
#if EIGEN_VERSION_AT_LEAST(3, 5, 0) // Eigen 5 reports itself as 3.5 in these macros
		// Eigen 5's BiCGSTAB compares the tolerance to the absolute residual norm (see EigenBiCGSTAB)
		_solver.setTolerance(this->Tolerance * b.norm());
#else
		_solver.setTolerance(this->Tolerance);
#endif
		if (xEquals0)
			x = _solver.solve(b);
		else
			x = _solver.solveWithGuess(b, x);
		this->IterationCount = _solver.iterations();
		this->LastIterationResult.NormalizedResidualNorm = (b - *this->Matrix * x).norm() / b.norm();
	}

private:
	void EndSetup(const SparseMatrix& A)
	{
		this->SetupComputationalWork = this->Precond.SetupComputationalWork();
		_solver.preconditioner().Precond = &this->Precond;
		_solver.compute(A);
	}
};
