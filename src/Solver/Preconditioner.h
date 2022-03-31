#pragma once
#include "IterativeSolver.h"
using namespace std;

class Preconditioner
{
public:
	Preconditioner() {}

	virtual void Setup(const SparseMatrix& A) {}
	virtual void Setup(const DenseMatrix& A) {}

	virtual Vector Apply(const Vector& r) = 0;

	virtual MFlops SetupComputationalWork() { return 0; }

	virtual MFlops SolvingComputationalWork() = 0;
};

class IdentityPreconditioner : public Preconditioner
{
public:
	Vector Apply(const Vector& r) override
	{
		return r;
	}
	MFlops SolvingComputationalWork() override
	{
		return 0;
	}
};

class DenseBlockJacobiPreconditioner : public Preconditioner
{
private:
	int _blockSize;
	vector<Eigen::FullPivLU<DenseMatrix>> _invD;
public:
	DenseBlockJacobiPreconditioner(int blockSize) 
	{
		_blockSize = blockSize;
	}

	void Setup(const DenseMatrix& A) override
	{
		auto nb = A.rows() / _blockSize;
		_invD = vector<Eigen::FullPivLU<DenseMatrix>>(nb);

		NumberParallelLoop<EmptyResultChunk> parallelLoop(nb);
		parallelLoop.Execute([this, &A](BigNumber i, ParallelChunk<EmptyResultChunk>* chunk)
			{
				DenseMatrix Di = A.block(i * _blockSize, i * _blockSize, _blockSize, _blockSize);
				_invD[i].compute(Di);
			});
	}

	MFlops SetupComputationalWork() override
	{
		return _invD.size() * 2.0 / 3.0 * pow(_blockSize, 3) * 1e-6;
	}

	Vector Apply(const Vector& r) override
	{
		int nb = _invD.size();
		assert(r.rows() == nb * _blockSize);

		Vector res(r.rows());
		for (int i = 0; i < nb; i++)
			res.segment(i * _blockSize, _blockSize) = _invD[i].solve(r.segment(i * _blockSize, _blockSize));
		return res;
	}

	MFlops SolvingComputationalWork() override
	{
		Utils::FatalError("TODO: SolvingComputationalWork() not implemented");
		return 0;
	}
};

class SolverPreconditioner : public Preconditioner
{
protected:
	Solver* _solver = nullptr;
public:
	SolverPreconditioner() {}

	SolverPreconditioner(Solver* solver)
	{
		_solver = solver;
		IterativeSolver* iterSolver = static_cast<IterativeSolver*>(solver);
		if (iterSolver)
		{
			iterSolver->ComputeExactSolution = false;
			iterSolver->PrintIterationResults = false;
			iterSolver->StoppingCrit = StoppingCriteria::MaxIterations;
			iterSolver->MaxIterations = 1;
		}
	}

	/*IterativeSolver* GetSolver()
	{
		return _solver;
	}*/

	/*bool IsIdentity()
	{
		return _solver == nullptr;
	}*/

	bool CanOptimizeResidualComputation()
	{
		if (!_solver) return false;
		IterativeSolver* iterSolver = static_cast<IterativeSolver*>(_solver);
		if (iterSolver)
			return iterSolver->CanOptimizeResidualComputation();
		return false;
	}

	friend ostream& operator<<(ostream& os, const SolverPreconditioner& p)
	{
		if (p._solver)
			os << *(p._solver);
		else
			os << "none";
		return os;
	}

	void Setup(const SparseMatrix& A) override
	{
		if (_solver)
			_solver->Setup(A);
	}

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F)
	{
		if (!_solver)
			return;
		IterativeSolver* iterSolver = static_cast<IterativeSolver*>(_solver);
		if (iterSolver)
			iterSolver->Setup(A, A_T_T, A_T_F, A_F_F);
	}

	Vector Apply(const Vector& r) override
	{
		if (_solver)
			return _solver->Solve(r);
		return r;
	}

	pair<Vector, Vector> ApplyAndComputeAe(const Vector& r)
	{
		assert(_solver);
		IterativeSolver* iterSolver = static_cast<IterativeSolver*>(_solver);
		assert(iterSolver);

		pair<Vector, Vector> p;
		auto&[e, Ae] = p;

		e = Vector::Zero(r.rows());
		bool eEquals0 = true;
		iterSolver->Solve(r, e, eEquals0, false, true);
		Ae = std::move(iterSolver->Ax);

		return p;
	}
	
	MFlops SetupComputationalWork() override
	{
		if (_solver)
			return _solver->SetupComputationalWork;
		return 0;
	}

	MFlops SolvingComputationalWork() override
	{
		if (_solver)
			return _solver->SolvingComputationalWork;
		return 0;
	}
};