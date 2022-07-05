#pragma once
#include <ios>
#include <fstream>
#include <unsupported/Eigen/SparseExtra>
#include "Utils.h"
#include "../Solver/IterationResult.h"
using namespace std;

class ExportModule
{
private:
	string _outputDirectory;
	string _filePrefix;
	string _vectorValueSeparator;
	string _iterValueSeparator = ",";

public:
	ExportModule() :
		ExportModule("")
	{}

	ExportModule(string outputDirectory) :
		ExportModule(outputDirectory, "")
	{}

	ExportModule(string outputDirectory, string filePrefix) :
		ExportModule(outputDirectory, filePrefix, ",")
	{}

	ExportModule(string outputDirectory, string filePrefix, string vectorValueSeparator) :
		_outputDirectory(outputDirectory),
		_filePrefix(filePrefix),
		_vectorValueSeparator(vectorValueSeparator)
	{}

	string OutputDirectory() const
	{
		return _outputDirectory;
	}

	void AddFilePrefix(string prefix)
	{
		_filePrefix += prefix;
	}

	string GetFilePathPrefix() const
	{
		if (_outputDirectory.empty())
			Utils::FatalError("output directory undefined");
		return _outputDirectory + "/" + _filePrefix;
	}
	string GetFilePath(string suffix, string extension) const
	{
		return GetFilePath(suffix + extension);
	}
	string GetFilePath(string suffix_with_extension) const
	{
		if (_filePrefix.empty())
			return GetFilePathPrefix() + suffix_with_extension;
		else
			return GetFilePathPrefix() + "_" + suffix_with_extension;
	}
	string GetDatFilePath(string suffix) const
	{
		return GetFilePath(suffix, ".dat");
	}

	void ExportMatrix(const SparseMatrix& M, string suffix) const
	{
		string filePath = GetDatFilePath(suffix);
		Eigen::saveMarket(M, filePath);
		cout << "Matrix exported: " << filePath << endl;
	}

	void ExportMatrix(const DenseMatrix& M, string suffix) const
	{
		string filePath = GetDatFilePath(suffix);
		Eigen::saveMarket(M, filePath);
		cout << "Matrix exported: " << filePath << endl;
	}

	void ExportVector(const Vector& v, string suffix) const
	{
		string filePath = GetDatFilePath(suffix);
		Eigen::saveMarketVector(v, filePath);
		cout << "Vector exported: " << filePath << endl;
	}

	void CleanFile(string suffix, string extension)
	{
		CleanFile(suffix + extension);
	}
	void CleanFile(string suffix_with_extension)
	{
		string filePath = GetFilePath(suffix_with_extension);
		remove(filePath.c_str());
	}

	void ExportNewVectorValue(double value, string suffix) const
	{
		string filePath = GetDatFilePath(suffix);
		ofstream file(filePath, ios_base::app | ios_base::out);

		file << std::scientific << value << _vectorValueSeparator; // << endl;
		file.close();
	}

	void Export(const IterationResult& result, string suffix)
	{
		string filePath = GetFilePath(suffix, ".csv");
		ofstream file(filePath, ios_base::app | ios_base::out);
		file.precision(2);
		string sep = _iterValueSeparator;

		if (result.IterationNumber == 0)
		{
			// Label row
			file << "Iter" << sep;
			file << "NormalizedRes" << sep;
			if (result.L2Error != -1)
				file << "L2Error" << sep;
			if (result.BoundaryL2Norm != -1)
				file << "BoundaryL2Norm" << sep;
			file << "IterConvRate" << sep;
			file << "AsympConvRate" << sep;
			file << "MatVec" << sep;
			file << "CPUTime";
			file << endl;
		}

		file << result.IterationNumber << sep;
		file << "\t" << std::scientific << result.NormalizedResidualNorm << sep;
		if (result.L2Error != -1)
			file << "\t" << std::scientific << result.L2Error << sep;
		if (result.BoundaryL2Norm != -1)
			file << "\t" << std::scientific << result.BoundaryL2Norm << sep;
		
		if (result.IterationNumber == 0)
			file << "\t\t" << " " << sep;
		else
			file << "\t" << std::defaultfloat << result.IterationConvRate << sep;
		
		if (result.IterationNumber == 0)
			file << "\t\t" << " " << sep;
		else
			file << "\t" << std::defaultfloat << result.AsymptoticConvRate << sep;
		
		file << "\t" << result.NumberOfFineMatVec() << sep;

		if (result.IterationNumber == 0)
			file << "\t" << "0";
		else
			file << "\t" << std::fixed << std::setprecision(1) << result.SolvingCPUTime().InSeconds();

		file << endl;
		file.close();
	}
};