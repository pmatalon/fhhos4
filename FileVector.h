#pragma once
#include<stdio.h>
#include<string>
using namespace std;

class FileVector
{
private:
	FILE* _file;
public:
	FileVector(string filePath)
	{
		this->_file = fopen(filePath.c_str(), "w");
	}
	void Add(double value)
	{
		fprintf(this->_file, "%.17e \n", value);
	}
	~FileVector()
	{
		fclose(this->_file);
	}
};

