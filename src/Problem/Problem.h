#pragma once
#include <string>
#include <unsupported/Eigen/SparseExtra>
#include "../Mesh/Mesh.h"
#include "../Utils/Utils.h"
using namespace std;

template <int Dim>
class Problem
{
protected:
	string _outputDirectory;
	string _fileName;
public:
	Mesh<Dim>* _mesh;
	SparseMatrix A;
	Vector b;
	Vector SystemSolution;

	Problem(Mesh<Dim>* mesh, string outputDirectory)
	{
		this->_mesh = mesh;
		this->_outputDirectory = outputDirectory;
	}

	string GetFilePath(string suffix)
	{
		return this->_outputDirectory + "/" + this->_fileName + "_" + suffix + ".dat";
	}

	virtual void PrintPhysicalProblem() = 0;

	void ExportMatrix(const SparseMatrix& M, string suffix)
	{
		string filePath = GetFilePath(suffix);
		Eigen::saveMarket(M, filePath);
	}

	void ExportVector(const Vector& M, string suffix)
	{
		string filePath = GetFilePath(suffix);
		Eigen::saveMarketVector(M, filePath);
	}

	virtual void ExtractSolution()
	{
		this->ExtractSolution(this->SystemSolution);
	}

	virtual ~Problem()
	{	}

protected:
	void ExtractSolution(Vector solution)
	{
		this->ExtractSolution(solution, "");
	}

	void ExtractSolution(Vector solution, string suffix)
	{
		string solutionFilePath = GetFilePath("solution" + suffix);
		Eigen::saveMarketVector(solution, solutionFilePath);
		cout << "Solution exported to \t" << solutionFilePath << endl;
	}

public:
	virtual void Assemble(Action action) = 0;

	virtual double L2Error(DomFunction exactSolution) = 0;

protected:
	virtual double L2Error(FunctionalBasis<Dim>* basis, const Vector& solution, DomFunction exactSolution)
	{
		struct ChunkResult
		{
			double absoluteError = 0;
			double normExactSolution = 0;
		};

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(_mesh->Elements);
		parallelLoop.Execute([basis, &solution, exactSolution](Element<Dim>* element, ParallelChunk<ChunkResult>* chunk)
			{
				auto approximate = basis->GetApproximateFunction(solution, element->Number * basis->NumberOfLocalFunctionsInElement(element));
				chunk->Results.absoluteError += element->L2ErrorPow2(approximate, exactSolution);
				chunk->Results.normExactSolution += element->Integral([exactSolution](DomPoint p) { return pow(exactSolution(p), 2); });
			});


		double absoluteError = 0;
		double normExactSolution = 0;

		parallelLoop.AggregateChunkResults([&absoluteError, &normExactSolution](ChunkResult chunkResult)
			{
				absoluteError += chunkResult.absoluteError;
				normExactSolution += chunkResult.normExactSolution;
			});

		absoluteError = sqrt(absoluteError);
		normExactSolution = sqrt(normExactSolution);
		return normExactSolution != 0 ? absoluteError / normExactSolution : absoluteError;
	}
};