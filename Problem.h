#pragma once
#include <string>
#include <Eigen/Sparse>
using namespace std;

class Problem
{
protected:
	string _solutionName;
public:
	Eigen::SparseMatrix<double> A;
	Eigen::VectorXd b;
	Eigen::VectorXd Solution;

	Problem(string solutionName)
	{
		this->_solutionName = solutionName;
	}

	void Solve()
	{
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.analyzePattern(this->A);
		solver.factorize(this->A);
		this->Solution = solver.solve(this->b);
	}

	virtual ~Problem()
	{	}
};