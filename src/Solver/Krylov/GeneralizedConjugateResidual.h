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
private:
	// The vectors of the iterations, allocated at the setup (AllocateIterationVectors()): the solve allocates nothing (see
	// Level::CycleVectors). Z and C: those of a cycle (Restart each).
	vector<Vector> _Z, _C;
	Vector _cycleCorrection;
public:
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
		AllocateIterationVectors(A.rows());
	}

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
		IterativeSolver::Setup(A);
		this->Precond.Setup(A, A_T_T, A_T_F, A_F_T, A_F_F, cellInterpOfOne, faceInterpOfOne);
		this->SetupComputationalWork = this->Precond.SetupComputationalWork();
		AllocateIterationVectors(A.rows());
	}

	void Setup(const DistributedMatrix& A) override
	{
		SetDistributedMatrix(A);
		this->Precond.Setup(A);
		this->SetupComputationalWork = this->Precond.SetupComputationalWork();
		AllocateIterationVectors(A.Diag.rows());
	}

	void Setup(const DistributedMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, shared_ptr<const IndexMap> faceMap, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
		SetDistributedMatrix(A);
		this->Precond.Setup(A, A_T_T, A_T_F, A_F_T, faceMap, cellInterpOfOne, faceInterpOfOne);
		this->SetupComputationalWork = this->Precond.SetupComputationalWork();
		AllocateIterationVectors(A.Diag.rows());
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
			r = ResidualOf(b, x);                                     result.AddWorkInFlops(Cost::DAXPY(A));
			result.SetResidualNorm(Parallel::Norm(r));                result.AddWorkInFlops(Cost::Norm(r));
		}

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			Utils::Log() << result << endl;

		// The current cycle: the preconditioned residuals z_j, the orthonormalized c_j, and Gamma (upper triangular) and a
		// such that A Z = C Gamma and r = r_0 - C a. x is the iterate at the start of the cycle.
		vector<Vector>& Z = _Z;
		vector<Vector>& C = _C;
		DenseMatrix Gamma = DenseMatrix::Zero(Restart, Restart);
		Vector a = Vector::Zero(Restart);
		int j = 0; // iterations done in the current cycle

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			// Apply the preconditioner: z_j = Prec(r), c_j = A z_j
			if (Precond.CanOptimizeResidualComputation())
			{
				Precond.ApplyAndComputeAe(r, Z[j], C[j]);             result.AddWorkInMFlops(Precond.SolvingComputationalWork());
			}
			else
			{
				Precond.Apply(r, Z[j]);                               result.AddWorkInMFlops(Precond.SolvingComputationalWork());
				if (this->Distributed)                                // C[j] = MatVec(Z[j]) in C[j]
					this->Distributed->Multiply(Z[j], nullptr, C[j]);
				else
					C[j].noalias() = A * Z[j];
				                                                      result.AddWorkInFlops(Cost::MatVec(A));
			}

			// Orthonormalization of c_j against the previous c_i (modified Gram-Schmidt). The vector operations are fused into
			// passes over the blocks of rows of the vectors, in parallel, each block updated then reduced while in cache
			// (Parallel::SumsOverRowBlocks()): the same operations on each coefficient as one operation after the other, and
			// the dot products of Parallel::Dot(). The pass i subtracts the projection on c_{i-1}, then computes c_i.c_j
			// (i < j) or |c_j|^2 (i = j).
			const Communicator& comm = Communicator::Current();
			BigNumber n = r.rows();
			Vector& c_j = C[j];
			for (int i = 0; i <= j; i++)
			{
				double sum = comm.Sum(Parallel::SumsOverRowBlocks<1>(n, [&](BigNumber first, BigNumber m)
				{
					if (i > 0)
						c_j.segment(first, m) -= Gamma(i - 1, j) * C[i - 1].segment(first, m);
					return array<double, 1>{ i < j ? Parallel::DotOfRows(C[i], c_j, first, m) : Parallel::SquaredNormOfRows(c_j, first, m) };
				})[0]);
				if (i < j)
				{
					Gamma(i, j) = sum;                                    result.AddWorkInFlops(Cost::Dot(c_j) + Cost::VectorDAXPY(c_j));
				}
				else
				{
					Gamma(j, j) = std::sqrt(sum);                         result.AddWorkInFlops(Cost::Norm(c_j));
				}
			}
			assert(Gamma(j, j) > 0);

			// Minimization of the residual norm along c_j: c_j normalized and a_j = c_j.r in one pass, then r and its norm
			a[j] = comm.Sum(Parallel::SumsOverRowBlocks<1>(n, [&](BigNumber first, BigNumber m)
			{
				c_j.segment(first, m) /= Gamma(j, j);
				return array<double, 1>{ Parallel::DotOfRows(c_j, r, first, m) };
			})[0]);                                                   result.AddWorkInFlops(Cost::ConstantVec(c_j) + Cost::Dot(r));
			double r_squaredNorm = comm.Sum(Parallel::SumsOverRowBlocks<1>(n, [&](BigNumber first, BigNumber m)
			{
				r.segment(first, m) -= a[j] * c_j.segment(first, m);
				return array<double, 1>{ Parallel::SquaredNormOfRows(r, first, m) };
			})[0]);                                                   result.AddWorkInFlops(Cost::VectorDAXPY(r));
			j++;

			//---------- End of iteration ----------//
			this->IterationCount++;

			result.SetResidualNorm(std::sqrt(r_squaredNorm));         result.AddWorkInFlops(Cost::Norm(r));

			if (j == Restart || StoppingCriteriaMet(result))
			{
				// End of the cycle: the iterate, then the restart from it
				CycleCorrection(Z, Gamma, a, j, _cycleCorrection, result);
				Parallel::ForEachRowBlock(x.rows(), [&](BigNumber first, BigNumber m) { x.segment(first, m) += _cycleCorrection.segment(first, m); });
				j = 0;

				// Residual replacement (the loop condition then evaluates the stopping criteria on the true residual)
				if (this->ResidualReplacement)
				{
					r = ResidualOf(b, x);                             result.AddWorkInFlops(Cost::DAXPY(A));
					result.ReplaceResidualNorm(Parallel::Norm(r));    result.AddWorkInFlops(Cost::Norm(r));
				}
				result.SetX(x);
			}
			else if (result.NeedsX()) // error against the exact solution, callback: the iterate at each iteration
			{
				CycleCorrection(Z, Gamma, a, j, _cycleCorrection, result);
				result.SetX(x + _cycleCorrection);
			}
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
	// Z Gamma^{-1} a, from the first j iterations of the cycle, in correction: the correction of the iterate of the start of
	// the cycle
	void CycleCorrection(const vector<Vector>& Z, const DenseMatrix& Gamma, const Vector& a, int j, Vector& correction, IterationResult& result) const
	{
		Vector y = Gamma.topLeftCorner(j, j).triangularView<Eigen::Upper>().solve(a.head(j));
		correction.resize(Z[0].rows());
		Parallel::ForEachRowBlock(Z[0].rows(), [&](BigNumber first, BigNumber m) // one pass (the same operations on each coefficient)
		{
			correction.segment(first, m) = y[0] * Z[0].segment(first, m);
			for (int i = 1; i < j; i++)
				correction.segment(first, m) += y[i] * Z[i].segment(first, m);
		});                                                           result.AddWorkInFlops(Cost::ConstantVec(Z[0]));
		for (int i = 1; i < j; i++)
		{
			                                                          result.AddWorkInFlops(Cost::VectorDAXPY(Z[i]));
		}
	}

	// r (Residual), the Restart vectors z_j and c_j of a cycle, and its correction: n rows
	void AllocateIterationVectors(BigNumber n)
	{
		this->Residual.resize(n);
		_Z.assign(Restart, Vector(n));
		_C.assign(Restart, Vector(n));
		_cycleCorrection.resize(n);
	}
};
