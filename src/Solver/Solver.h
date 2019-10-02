#pragma once
#include "../Utils/Utils.h"
using namespace std;

class Solver
{
public:
	BigNumber SetupComputationalWork = 0;
	BigNumber SolvingComputationalWork = 0;

	Solver() { }

	virtual void Serialize(ostream& os) const = 0;

	friend ostream& operator<<(ostream& os, const Solver& s)
	{
		s.Serialize(os);
		return os;
	}

	virtual void Setup(const SparseMatrix& A) = 0;

	virtual Eigen::VectorXd Solve(const Eigen::VectorXd& b) = 0;

	virtual ~Solver() {}
};