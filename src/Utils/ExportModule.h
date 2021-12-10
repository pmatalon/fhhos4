#pragma once
#include <string>
#include <unsupported/Eigen/SparseExtra>
#include "Utils.h"
using namespace std;

class ExportModule
{
private:
	string _outputDirectory;
	string _filePrefix;

public:
	ExportModule() :
		ExportModule("")
	{}

	ExportModule(string outputDirectory) :
		ExportModule(outputDirectory, "")
	{}

	ExportModule(string outputDirectory, string filePrefix) :
		_outputDirectory(outputDirectory),
		_filePrefix(filePrefix)
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
		if (_filePrefix.empty())
			return GetFilePathPrefix() + suffix + extension;
		else
			return GetFilePathPrefix() + "_" + suffix + extension;
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

	void ExportVector(const Vector& v, string suffix) const
	{
		string filePath = GetDatFilePath(suffix);
		Eigen::saveMarketVector(v, filePath);
		cout << "Vector exported: " << filePath << endl;
	}


};