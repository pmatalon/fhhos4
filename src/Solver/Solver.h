#pragma once
#include "../Utils/Utils.h"
#include "../Utils/Cost.h"
using namespace std;

class Solver
{
public:
	const SparseMatrix* Matrix = nullptr;

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
	}

	virtual Vector Solve(const Vector& b) = 0;

	virtual ~Solver() {}
};