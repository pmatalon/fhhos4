#pragma once
#include "../Utils/Utils.h"
using namespace std;

class Solver
{
public:
	size_t SetupComputationalWork = 0;
	size_t SolvingComputationalWork = 0;

	Solver() { }

	virtual void Serialize(ostream& os) const = 0;

	friend ostream& operator<<(ostream& os, const Solver& s)
	{
		s.Serialize(os);
		return os;
	}

	virtual void Setup(const SparseMatrix& A) = 0;

	virtual Vector Solve(const Vector& b) = 0;

	virtual ~Solver() {}
};