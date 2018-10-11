#pragma once
#include<stdio.h>
#include<string>
using namespace std;

class FileMatrix
{
private:
	int _nRows;
	int _nCols;
	FILE* _file;
public:
	FileMatrix(int nRows, int nCols, string filePath);
	void Add(int i, int j, double value);
	~FileMatrix();
};

FileMatrix::FileMatrix(int nRows, int nCols, string filePath)
{
	this->_nRows = nRows;
	this->_nCols = nCols;
	this->_file = fopen(filePath.c_str(), "w");
}

void FileMatrix::Add(int i, int j, double value)
{
	if (value == 0)
		return;
	fprintf(this->_file, "%llu %llu %.17g \n", i, j, value);
}

FileMatrix::~FileMatrix()
{
	this->Add(this->_nRows, this->_nCols, 0);
	fclose(this->_file);
}

