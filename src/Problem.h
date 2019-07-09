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