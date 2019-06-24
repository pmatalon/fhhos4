#pragma once
#include <string>
#include <Eigen/Sparse>
#include "Utils/agmg.h"
using namespace std;

#ifdef TPL_ENABLE_AGMG
extern "C" void dagmg_(int*,double*,int*,int*,double*,double*,int*,int*,int*,int*,double*);
#endif

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

	virtual void Solve()
	{
		cout << "Solving..." << endl;
#ifdef TPL_ENABLE_AGMG
		int ijob=0,nrest=1;
		int nun = this->A.cols();
		int nnz = this->A.nonZeros();

		/*      maximal number of iterations */
		int iter=100;
		/*      tolerance on relative residual norm */
		double tol=1.e-10;
		/*      unit number for output messages: 6 => standard output */
		int iprint=6;

		/*      generate the matrix in required format (CSR) */
		/*        first allocate the vectors with correct size */
		double *a,*f,*x;
		int *ja,*ia;
		ia = new int[nun+1];
		ja = new int[nnz];
		a = new double[nnz];
		f = new double[nun];
		x = new double[nun];

		/*        next set entries */
		for(int i = 0; i < nnz; i++)
		{
			a[i] = this->A.valuePtr()[i];
			ja[i] = this->A.innerIndexPtr()[i] + 1;
		}
		for(int i = 0; i < nun; i++)
		{
			ia[i] = this->A.outerIndexPtr()[i] + 1;
			f[i] = this->b[i];
		}
		ia[nun] = this->A.outerIndexPtr()[nun] + 1; // last value

		/*      call agmg
		argument 5 (ijob)  is 0 because we want a complete solve
		argument 7 (nrest) is 1 because we want to use flexible CG
		                 (the matrix is symmetric positive definite) */
		dagmg_(&nun,a,ja,ia,f,x,&ijob,&iprint,&nrest,&iter,&tol);

		// copy into solution
		this->Solution = Eigen::VectorXd(nun);
		for(int i = 0; i < nun; i++)
			this->Solution[i] = x[i];

		delete ia, ja, a, f, x;
#else
		/*Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.analyzePattern(this->A);
		solver.factorize(this->A);*/

		Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double>> solver;
		solver.setTolerance(1.e-10);
		solver.compute(this->A);

		this->Solution = solver.solve(this->b);
#endif
	}

	string GetFilePath(string suffix)
	{
		return this->_outputDirectory + "/" + this->_fileName + "_" + suffix + ".dat";
	}

	void ExportMatrix(const Eigen::SparseMatrix<double>& M, string suffix)
	{
		string filePath = GetFilePath(suffix);
		Eigen::saveMarket(M, filePath);
	}

	void ExportVector(const Eigen::VectorXd& M, string suffix)
	{
		string filePath = GetFilePath(suffix);
		Eigen::saveMarketVector(M, filePath);
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
		this->ExtractSolution(solution, "");
	}

	void ExtractSolution(Eigen::VectorXd solution, string suffix)
	{
		string solutionFilePath = GetFilePath("solution" + suffix);
		Eigen::saveMarketVector(solution, solutionFilePath);
		cout << "Solution exported to \t" << solutionFilePath << endl;
	}
};