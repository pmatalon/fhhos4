#pragma once
#include <string>
#include <Eigen/Sparse>
using namespace std;

class Problem
{
protected:
	string _solutionName;
	string _outputDirectory;
	string _fileName;
public:
	Eigen::SparseMatrix<double> A;
	Eigen::VectorXd b;
	Eigen::VectorXd Solution;

	Problem(string solutionName, string outputDirectory)
	{
		this->_solutionName = solutionName;
		this->_outputDirectory = outputDirectory;
	}

	void Solve()
	{
		cout << "Solving..." << endl;
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.analyzePattern(this->A);
		solver.factorize(this->A);

		/*Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double>> solver;
		solver.setTolerance(1.e-10);
		solver.compute(this->A);*/

		this->Solution = solver.solve(this->b);
	}

	virtual void ExtractSolution()
	{
		this->ExtractSolution(this->Solution);
	}

	virtual ~Problem()
	{	}

protected:
	void ExtractSolution(Eigen::VectorXd solution)
	{
		string solutionFilePath = this->_outputDirectory + "/" + this->_fileName + "_solution.dat";
		Eigen::saveMarketVector(solution, solutionFilePath);
		cout << "Solution exported to \t" << solutionFilePath << endl;
	}
};