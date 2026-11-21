#pragma once
#include "../IterativeSolver.h"
#include "../Preconditioner.h"
#include <list>
#include <memory>
using namespace std;

class FlexibleConjugateGradient : public IterativeSolver
{
private:
	struct Direction
	{
		Vector* d;
		Vector* Ad;
		double d_dot_Ad;
		double Coefficient = 0; // in the orthogonalization of the current iteration: z.Ad / d.Ad
	};
	// The vectors of the iterations, allocated at the setup (AllocateIterationVectors()): the solve allocates nothing (a
	// large vector would be fresh memory at each allocation, see Level::CycleVectors). The spare vectors: those of the
	// directions not in use, for the next ones.
	Vector _z;
	vector<unique_ptr<Vector>> _spareVectors;
public:
	SolverPreconditioner Precond;
	int Truncation; // max number of previous search direction to use for A-orthogonalization
	bool Restart = false;
	// Residual replacement: every ResidualRecomputation iterations, and when the stopping criteria are met, the residual
	// updated by the recurrence (r -= alpha*Ad) is replaced by the true residual b - Ax, on which the stopping criteria
	// are then evaluated (the iterations continue if they are no longer met). In finite precision, the recurrence drifts
	// from the true residual near the attainable accuracy (~eps |A| |x|), below which the updated residual still
	// decreases (as hypre's PCG options RecomputeResidualP and RecomputeResidual; van der Vorst and Ye, SIAM J. Sci.
	// Comput. 22(3), 2000). 0: never (the program: the iteration counts of the papers).
	int ResidualRecomputation = 0;

	FlexibleConjugateGradient(int truncation = 1, bool restart = false)
	{
		Truncation = truncation;
		Restart = restart;
	}

	virtual void Serialize(ostream& os) const override
	{
		os << "Flexible Conjugate Gradient, preconditioner: " << Precond;
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
		if (A_F_T.rows() > 0)
			Utils::FatalError("The flexible conjugate gradient requires a symmetric matrix: use GCR for a non-symmetric one.");
		IterativeSolver::Setup(A);
		this->Precond.Setup(A, A_T_T, A_T_F, A_F_T, A_F_F, cellInterpOfOne, faceInterpOfOne);
		this->SetupComputationalWork = this->Precond.SetupComputationalWork();
		AllocateIterationVectors(A.rows());
	}

	// On several processes: A is read in full (both triangles)
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

public:
	void Solve(const Vector& b, Vector& x, bool xEquals0, bool computeResidual, bool computeAx) override
	{
		const SparseMatrix& A = *this->Matrix;

		// Compared to CG, the work of 1 iteration of FCG is 1 more dot product by previous direction saved.
		// FCG(1) costs 1 more dot product by iteration than CG, i.e. 1 matrix-vector multiplication + 3 dot products.
		this->SolvingComputationalWork = 0; 

		if (this->ComputeExactSolution)
			this->ExactSolution = this->_directSolver.solve(b);

		IterationResult result = CreateFirstIterationResult(b, x);

		// Init restart parameters
		list<Direction> previousDirections;

		Vector& r = this->Residual;
		if (xEquals0)
		{
			r = b;
			result.SetResidualAsB();
		}
		else
		{
			r = ResidualOfSymmetric(b, x);                            result.AddWorkInFlops(Cost::DAXPY(A));
			result.SetResidualNorm(Parallel::Norm(r));                result.AddWorkInFlops(Cost::Norm(r));
		}

		this->IterationCount = 0;

		if (this->PrintIterationResults)
			Utils::Log() << result << endl;

		while (!StoppingCriteriaReached(result))
		{
			result = IterationResult(result);

			// Apply the preconditioner (A z: the A d of the new direction, before its orthogonalization)
			Vector& z = _z;
			Vector* Ad = nullptr;
			if (Precond.CanOptimizeResidualComputation())
			{
				Ad = SpareVector();
				Precond.ApplyAndComputeAe(r, z, *Ad);              result.AddWorkInMFlops(Precond.SolvingComputationalWork());
			}
			else
			{
				Precond.Apply(r, z);                               result.AddWorkInMFlops(Precond.SolvingComputationalWork());
			}

			// The final direction of research, d, is found by A-orthogonalizing z with the previous directions of research,
			// according to the truncation / restart method chosen: d = z - sum_k (z.Ad_k / d_k.Ad_k) d_k, and A d from A z
			// alike. The vector operations are fused into passes over the blocks of rows of the vectors, in parallel, each
			// block updated then reduced while in cache (Parallel::SumsOverRowBlocks()): the same operations on each
			// coefficient as one operation after the other, and the dot products of Parallel::Dot().
			const Communicator& comm = Communicator::Current();
			BigNumber n = r.rows();
			assert(x.rows() == n);
			bool AdOfPreconditioner = Ad != nullptr;
			for (Direction& directionk : previousDirections)
			{
				double z_dot_Adk = Parallel::Dot(z, *directionk.Ad);      result.AddWorkInFlops(Cost::Dot(z));
				directionk.Coefficient = z_dot_Adk / directionk.d_dot_Ad; result.AddWorkInFlops(Cost::VectorDAXPY(z));
				if (AdOfPreconditioner)
				{
					                                                      result.AddWorkInFlops(Cost::VectorDAXPY(z));
				}
			}
			Vector* d = SpareVector();
			d->resize(n);
			auto orthogonalize = [&](BigNumber i, BigNumber m)
			{
				d->segment(i, m) = z.segment(i, m);
				for (const Direction& directionk : previousDirections)
				{
					d->segment(i, m) -= directionk.Coefficient * directionk.d->segment(i, m);
					if (AdOfPreconditioner)
						Ad->segment(i, m) -= directionk.Coefficient * directionk.Ad->segment(i, m);
				}
			};
			double d_dot_Ad;
			double r_dot_d;
			if (AdOfPreconditioner)
			{
				array<double, 2> dots = Parallel::SumsOverRowBlocks<2>(n, [&](BigNumber i, BigNumber m)
				{
					orthogonalize(i, m);
					return array<double, 2>{ Parallel::DotOfRows(*d, *Ad, i, m), Parallel::DotOfRows(r, *d, i, m) };
				});
				d_dot_Ad = comm.Sum(dots[0]);
				r_dot_d = comm.Sum(dots[1]);
			}
			else
			{
				r_dot_d = comm.Sum(Parallel::SumsOverRowBlocks<1>(n, [&](BigNumber i, BigNumber m)
				{
					orthogonalize(i, m);
					return array<double, 1>{ Parallel::DotOfRows(r, *d, i, m) };
				})[0]);
				Ad = SpareVector();
				if (this->Distributed)
					this->Distributed->Multiply(*d, nullptr, *Ad);
				else
					Ad->noalias() = A.selfadjointView<Eigen::Lower>() * (*d);
				                                                      result.AddWorkInFlops(Cost::MatVec(A));
				d_dot_Ad = Parallel::Dot(*d, *Ad);
			}
			                                                          result.AddWorkInFlops(Cost::Dot(*d));
			assert(d_dot_Ad != 0);

			// Step length in the direction of research
			double alpha = r_dot_d / d_dot_Ad;                        result.AddWorkInFlops(Cost::Dot(r));

			// Moving from the current solution to the next by taking the step in the direction of research, the new
			// residual and its norm: one pass
			double r_squaredNorm = comm.Sum(Parallel::SumsOverRowBlocks<1>(n, [&](BigNumber i, BigNumber m)
			{
				x.segment(i, m) += alpha * d->segment(i, m);
				r.segment(i, m) -= alpha * Ad->segment(i, m);
				return array<double, 1>{ Parallel::SquaredNormOfRows(r, i, m) };
			})[0]);
			                                                          result.AddWorkInFlops(Cost::VectorDAXPY(x));
			                                                          result.AddWorkInFlops(Cost::VectorDAXPY(r));

			// Restart?
			if (previousDirections.size() == this->Truncation)
			{
				if (this->Restart)
					Clear(previousDirections);
				else
					DeleteOldest(previousDirections);
			}

			// The direction of research is saved 
			// (even after a restart, as advised by Notay)
			Direction direction = { d, Ad, d_dot_Ad };
			previousDirections.push_back(direction);

			//---------- End of iteration ----------//
			this->IterationCount++;

			result.SetX(x);
			result.SetResidualNorm(std::sqrt(r_squaredNorm));         result.AddWorkInFlops(Cost::Norm(r));

			// Residual replacement (the loop condition then evaluates the stopping criteria on the true residual)
			if (this->ResidualRecomputation > 0 && (this->IterationCount % this->ResidualRecomputation == 0 || StoppingCriteriaMet(result)))
			{
				r = ResidualOfSymmetric(b, x);                        result.AddWorkInFlops(Cost::DAXPY(A));
				result.ReplaceResidualNorm(Parallel::Norm(r));        result.AddWorkInFlops(Cost::Norm(r));
			}

			if (this->PrintIterationResults)
				Utils::Log() << result << endl;
		}

		Clear(previousDirections);

		if (this->PrintIterationResults)
			Utils::Log() << endl;

		result.CopyInfoInto(this->LastIterationResult);
		this->SolvingComputationalWork = result.SolvingComputationalWork();
	}

private:
	// b - A x, A read from its lower triangular part (symmetric), or in full on several processes
	Vector ResidualOfSymmetric(const Vector& b, const Vector& x) const
	{
		if (this->Distributed)
			return b - *this->Distributed * x;
		return b - this->Matrix->selfadjointView<Eigen::Lower>() * x;
	}

	// The vectors of the iterations, of n rows: r (Residual), z, and those of the Truncation + 1 directions in use at most
	// (the previous ones and the new one)
	void AllocateIterationVectors(BigNumber n)
	{
		this->Residual.resize(n);
		_z.resize(n);
		_spareVectors.clear();
		for (int i = 0; i < 2 * (Truncation + 1); i++)
			_spareVectors.push_back(make_unique<Vector>(n));
	}

	// A vector for a new direction: a spare one (allocated at the setup)
	Vector* SpareVector()
	{
		if (_spareVectors.empty())
			return new Vector();
		Vector* v = _spareVectors.back().release();
		_spareVectors.pop_back();
		return v;
	}

	void Clear(list<Direction>& previousDirections)
	{
		for (Direction& directionk : previousDirections)
		{
			_spareVectors.emplace_back(directionk.d);
			_spareVectors.emplace_back(directionk.Ad);
		}
		previousDirections.clear();
	}

	void DeleteOldest(list<Direction>& previousDirections)
	{
		Direction& oldestDirection = previousDirections.front();
		_spareVectors.emplace_back(oldestDirection.d);
		_spareVectors.emplace_back(oldestDirection.Ad);
		previousDirections.pop_front();
	}
};