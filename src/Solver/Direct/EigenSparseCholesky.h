#pragma once
#include "../Solver.h"
#include "GatheredOnRoot.h"
#include "MUMPS.h"
#include <Eigen/SparseCholesky>
using namespace std;

class EigenSparseCholesky : public Solver
{
private:
	// Eigen's sparse Cholesky requires a column-major matrix type.
	// Memory: SimplicialLDLT only keeps the factors (it makes a temporary permuted copy of the matrix in compute()).
	Eigen::SimplicialLDLT<ColMajorSparseMatrix> _solver;
	bool _gathered; // on several processes: the matrix gathered on the process 0 (the reproducible modes), or MUMPS's
	unique_ptr<MUMPSSolver> _distributed;
	unique_ptr<GatheredOnRoot> _gatheredOnRoot;

public:
	explicit EigenSparseCholesky(bool gathered = false) : Solver(), _gathered(gathered) {}

	void Serialize(ostream& os) const override
	{
		os << "Cholesky factorization (Eigen library; on several processes: " << (_gathered ? "Eigen's on the process 0" : "MUMPS's LDL^T") << ")";
	}

	void Setup(const SparseMatrix& A) override
	{
		Solver::Setup(A);
		_distributed.reset();
		_gatheredOnRoot.reset();
		// Copy: A is row-major, SimplicialLDLT needs column-major. Freed at the end of this function.
		// Since A is symmetric, we copy its transpose: it is the same matrix, and the transpose of
		// a row-major matrix is column-major, so the entries are copied in order.
		ColMajorSparseMatrix colMajorA = A.transpose();
		_solver.compute(colMajorA);
		this->SetupComputationalWork = Cost::CholeskyFactorization(A)*1e-6;
		Eigen::ComputationInfo info = _solver.info();
		if (info != Eigen::ComputationInfo::Success)
		{
			//cout << "----------------- A -------------------" << A << endl;
			Utils::FatalError("SimplicialLDLT failed to execute with the code " + to_string(info) + ".");
		}
	}

	// On several processes (a coarse solver): MUMPS's LDL^T, distributed, or the matrix gathered on the process 0 and
	// factorized there (the reproducible modes). Collective.
	void Setup(const DistributedMatrix& A) override
	{
		if (_gathered)
			_gatheredOnRoot = make_unique<GatheredOnRoot>(A, [this](const SparseMatrix& M) { Setup(M); });
		else
		{
			_distributed = make_unique<MUMPSSolver>(true);
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
		this->SolvingComputationalWork = Cost::CholeskySolve(_solver.matrixL())*1e-6;
		return x;
	}
};