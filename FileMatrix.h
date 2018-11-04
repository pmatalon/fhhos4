
#include<stdio.h>
#include<string>
using namespace std;

#pragma once
class FileMatrix
{
private:
	BigNumber _nRows;
	BigNumber _nCols;
	FILE* _file;
public:
	FileMatrix(BigNumber nRows, BigNumber nCols, string filePath)
	{
		this->_nRows = nRows;
		this->_nCols = nCols;
		this->_file = fopen(filePath.c_str(), "w");
	}
	void Add(BigNumber i, BigNumber j, double value)
	{
		if (value == 0)
			return;
		fprintf(this->_file, "%llu %llu %.17g \n", i, j, value);
	}
	~FileMatrix()
	{
		this->Add(this->_nRows, this->_nCols, 0);
		fclose(this->_file);
	}
};
