#pragma once
#include<stdio.h>
#include<string>
using namespace std;

class FileVector
{
private:
	FILE* _file;
public:
	FileVector(string filePath);
	void Add(double value);
	~FileVector();
};

FileVector::FileVector(string filePath)
{
	this->_file = fopen(filePath.c_str(), "w");
}

void FileVector::Add(double value)
{
	if (value == 0)
		return;
	fprintf(this->_file, "%.17g \n", value);
}

FileVector::~FileVector()
{
	fclose(this->_file);
}

