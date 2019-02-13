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
		cout << "Solving..." << endl;
		/*Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.analyzePattern(this->A);
		solver.factorize(this->A);*/

		Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double>> solver;
		solver.setTolerance(1.e-10);
		solver.compute(this->A);

		this->Solution = solver.solve(this->b);
	}

	virtual ~Problem()
	{	}
};