#pragma once
#include "../IterativeSolver.h"
#include "../Preconditioner.h"
using namespace std;

// Preconditioned GCR (generalized conjugate residual), restarted: GCR(Restart), for non-symmetric matrices. Minimizes
// the residual norm over the space of the preconditioned residuals, and allows a variable (non-linear) preconditioner,
// as the K-cycle. Economical version of Y. Notay, "An aggregation-based algebraic multigrid method", ETNA 37 (2010),
// Algorithm 3.1: only the vectors c_j = A z_j are orthonormalized (modified Gram-Schmidt), and the iterate
// x = x_0 + Z Gamma^{-1} a is formed once per cycle (restart or end), which costs as much as FGMRES per iteration,
// without its matrix-vector product at the restart. Reads the whole matrix (no symmetry assumed).
class GeneralizedConjugateResidual : public IterativeSolver
{
public:
	SolverPreconditioner Precond;
	int Restart; // maximal number of iterations of a cycle
	// Residual replacement: at the end of each cycle (restart, or stopping criteria met), the residual updated by the
	// recurrence (r -= alpha*c) is replaced by the true residual b - Ax, on which the stopping criteria are then evaluated
	// (the iterations continue if they are no longer met). See FlexibleConjugateGradient::ResidualRecomputation. Costs one
	// matrix-vector product per cycle. false: Notay's algorithm.
	bool ResidualReplacement = false;

	GeneralizedConjugateResidual(int restart = 10)
	{
		Restart = restart;
	}

	virtual void Serialize(ostream& os) const override
	{
		os << "GCR(" << Restart << "), preconditioner: " << Precond;
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		this->Precond.Setup(A);
		this->SetupComputationalWork = this->Precond.SetupComputationalWork();
	}

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
		IterativeSolver::Setup(A);
		this->Precond.Setup(A, A_T_T, A_T_F, A_F_T, A_F_F, cellInterpOfOne, faceInterpOfOne);
		this->SetupComputationalWork = this->Precond.SetupComputationalWork();
	}

	void Solve(const Vector& b, Vector& x, bool xEquals0, bool computeResidual, bool computeAx) override
	{
		const SparseMatrix& A = *this->Matrix;

		this->SolvingComputationalWork = 0;

		if (this->ComputeExactSolution)
			this->ExactSolution = this->_directSolver.solve(b);

		IterationResult result = CreateFirstIterationResult(b, x);

		Vector& r = this->Residual;
		if (xEquals0)
		{
			r = b;
			result.SetResidualAsB();
		}
		else
		{
			r = b - A * x;                                            result.AddWorkInFlops(Cost::DAXPY(A));
			result.SetResidualNorm(r.norm());                         result.AddWorkInFlops(Cost::Norm(r));
		}

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			Utils::Log() << result << endl;

		// The current cycle: the preconditioned residuals z_j, the orthonormalized c_j, and Gamma (upper triangular) and a
		// such that A Z = C Gamma and r = r_0 - C a. x is the iterate at the start of the cycle.
		vector<Vector> Z(Restart), C(Restart);
		DenseMatrix Gamma = DenseMatrix::Zero(Restart, Restart);
		Vector a = Vector::Zero(Restart);
		int j = 0; // iterations done in the current cycle

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			// Apply the preconditioner: z_j = Prec(r), c_j = A z_j
			if (Precond.CanOptimizeResidualComputation())
			{
				std::tie(Z[j], C[j]) = Precond.ApplyAndComputeAe(r);  result.AddWorkInMFlops(Precond.SolvingComputationalWork());
			}
			else
			{
				Z[j] = Precond.Apply(r);                              result.AddWorkInMFlops(Precond.SolvingComputationalWork());
				C[j] = A * Z[j];                                      result.AddWorkInFlops(Cost::MatVec(A));
			}

			// Orthonormalization of c_j against the previous c_i (modified Gram-Schmidt)
			for (int i = 0; i < j; i++)
			{
				Gamma(i, j) = C[i].dot(C[j]);                         result.AddWorkInFlops(Cost::Dot(C[j]));
				C[j] -= Gamma(i, j) * C[i];                           result.AddWorkInFlops(Cost::VectorDAXPY(C[j]));
			}
			Gamma(j, j) = C[j].norm();                                result.AddWorkInFlops(Cost::Norm(C[j]));
			assert(Gamma(j, j) > 0);
			C[j] /= Gamma(j, j);                                      result.AddWorkInFlops(Cost::ConstantVec(C[j]));

			// Minimization of the residual norm along c_j
			a[j] = C[j].dot(r);                                       result.AddWorkInFlops(Cost::Dot(r));
			r -= a[j] * C[j];                                         result.AddWorkInFlops(Cost::VectorDAXPY(r));
			j++;

			//---------- End of iteration ----------//
			this->IterationCount++;

			result.SetResidualNorm(r.norm());                         result.AddWorkInFlops(Cost::Norm(r));

			if (j == Restart || StoppingCriteriaMet(result))
			{
				// End of the cycle: the iterate, then the restart from it
				x += CycleCorrection(Z, Gamma, a, j, result);
				j = 0;

				// Residual replacement (the loop condition then evaluates the stopping criteria on the true residual)
				if (this->ResidualReplacement)
				{
					r = b - A * x;                                    result.AddWorkInFlops(Cost::DAXPY(A));
					result.ReplaceResidualNorm(r.norm());             result.AddWorkInFlops(Cost::Norm(r));
				}
				result.SetX(x);
			}
			else if (result.NeedsX()) // error against the exact solution, callback: the iterate at each iteration
				result.SetX(x + CycleCorrection(Z, Gamma, a, j, result));
			else
				result.StopTimers();

			if (this->PrintIterationResults)
				Utils::Log() << result << endl;
		}

		if (this->PrintIterationResults)
			Utils::Log() << endl;

		result.CopyInfoInto(this->LastIterationResult);
		this->SolvingComputationalWork = result.SolvingComputationalWork();
	}

private:
	// Z Gamma^{-1} a, from the first j iterations of the cycle: the correction of the iterate of the start of the cycle
	Vector CycleCorrection(const vector<Vector>& Z, const DenseMatrix& Gamma, const Vector& a, int j, IterationResult& result) const
	{
		Vector y = Gamma.topLeftCorner(j, j).triangularView<Eigen::Upper>().solve(a.head(j));
		Vector correction = y[0] * Z[0];                              result.AddWorkInFlops(Cost::ConstantVec(Z[0]));
		for (int i = 1; i < j; i++)
		{
			correction += y[i] * Z[i];                                result.AddWorkInFlops(Cost::VectorDAXPY(Z[i]));
		}
		return correction;
	}
};
