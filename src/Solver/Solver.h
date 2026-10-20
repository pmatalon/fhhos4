#pragma once
#include "../Utils/Utils.h"
#include "../Utils/Cost.h"
#include "../Utils/DistributedMatrix.h"
using namespace std;

class Solver
{
public:
	const SparseMatrix* Matrix = nullptr;
	const DenseMatrix* DenseMat = nullptr;
	// On several processes (MPI): the matrix of the setup, distributed by rows (Matrix is then its block Diag); null on
	// one process
	const DistributedMatrix* Distributed = nullptr;

	MFlops SetupComputationalWork = 0;
	MFlops SolvingComputationalWork = 0;

	Solver() { }

	virtual void Serialize(ostream& os) const = 0;

	friend ostream& operator<<(ostream& os, const Solver& s)
	{
		s.Serialize(os);
		return os;
	}

	virtual void Setup(const SparseMatrix& A)
	{
		this->Matrix = &A;
		this->Distributed = nullptr;
	}

	virtual void Setup(const DenseMatrix& A)
	{
		this->DenseMat = &A;
	}

	// Setup on several processes, for the solvers that run on them (U-AMG and what it uses)
	virtual void Setup(const DistributedMatrix& A)
	{
		Utils::FatalError("This solver does not run on several processes (MPI).");
	}

	virtual Vector Solve(const Vector& b) = 0;

	virtual Vector Solve()
	{
		Utils::FatalError("Solve() not implemented in this solver.");
		return Vector::Zero(0);
	}

	virtual ~Solver() {}
};