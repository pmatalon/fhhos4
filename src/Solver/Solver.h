#pragma once
#include <Eigen/Sparse>
using namespace std;

class Solver
{
public:
	Solver() { }

	virtual void Serialize(ostream& os) const = 0;

	friend ostream& operator<<(ostream& os, const Solver& s)
	{
		s.Serialize(os);
		return os;
	}

	virtual void Setup(const Eigen::SparseMatrix<double>& A) = 0;

	virtual Eigen::VectorXd Solve(const Eigen::VectorXd& b) = 0;

	virtual ~Solver() {}
};