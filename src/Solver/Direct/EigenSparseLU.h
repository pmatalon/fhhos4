#pragma once
#include "../Solver.h"
#include "GatheredOnRoot.h"
#include "MUMPS.h"
using namespace std;

class EigenSparseLU : public Solver
{
private:
	// Eigen's SparseLU requires a column-major matrix type. With a row-major one, it reads the row pointers
	// as column pointers and silently factorizes a wrong matrix when the sparsity pattern is not symmetric.
	// Memory: besides the L and U factors, SparseLU keeps its own (column-permuted) copy of the matrix.
	Eigen::SparseLU<ColMajorSparseMatrix> _solver;
	bool _gathered; // on several processes: the matrix gathered on the process 0 (the reproducible modes), or MUMPS's
	unique_ptr<MUMPSSolver> _distributed;
	unique_ptr<GatheredOnRoot> _gatheredOnRoot;

public:
	explicit EigenSparseLU(bool gathered = false) : Solver(), _gathered(gathered) {}

	void Serialize(ostream& os) const override
	{
		os << "LU factorization (Eigen library; on several processes: " << (_gathered ? "Eigen's on the process 0" : "MUMPS's") << ")";
	}

	void Setup(const SparseMatrix& A) override
	{
		Solver::Setup(A);
		_distributed.reset();
		_gatheredOnRoot.reset();
		//_solver.isSymmetric(true);
		// Copy: A is row-major, SparseLU needs column-major (see above). Freed at the end of this function.
		ColMajorSparseMatrix colMajorA = A;
		_solver.compute(colMajorA);
		this->SetupComputationalWork = Cost::LUFactorization(A)*1e-6;
		Eigen::ComputationInfo info = _solver.info();
		if (info != Eigen::ComputationInfo::Success)
		{
			//cout << "----------------- A -------------------" << A << endl;
			Utils::FatalError("SparseLU failed to execute with the code " + to_string(info) + ": " + _solver.lastErrorMessage());
		}
	}

	// On several processes (a coarse solver): MUMPS's LU, distributed, or the matrix gathered on the process 0 and
	// factorized there (the reproducible modes). Collective.
	void Setup(const DistributedMatrix& A) override
	{
		if (_gathered)
			_gatheredOnRoot = make_unique<GatheredOnRoot>(A, [this](const SparseMatrix& M) { Setup(M); });
		else
		{
			_distributed = make_unique<MUMPSSolver>(false);
			_distributed->Setup(A);
		}
		this->Matrix = &A.Diag;
		this->Distributed = &A;
	}

	Vector Solve(const Vector& b) override
	{
		if (_distributed)
			return _distributed->Solve(b);
		if (_gatheredOnRoot)
			return _gatheredOnRoot->SolveWith(b, [this](const Eigen::Map<const Vector>& all) { return Vector(_solver.solve(all)); });
		Vector x = _solver.solve(b);
		this->SolvingComputationalWork = 0;//Cost::LUSolve(_solver.matrixL()., _solver.matrixU()); TODO
		return x;
	}
};