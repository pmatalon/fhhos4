#pragma once
#include <fstream>
#include <sys/stat.h>
using namespace std;

#ifndef ROOT_PATH
#define ROOT_PATH "./"
#endif // !ROOT_PATH

class FileSystem
{
public:
	static string RootPath()
	{
		return ROOT_PATH;
	}

	static bool FileExists(string filename)
	{
		ifstream ifile(filename);
		return ifile.good();
	}

	static string FileName(const string& filePath)
	{
		return filePath.substr(filePath.find_last_of("/\\") + 1);
	}

	static bool HasExtension(const string& filePath)
	{
		return filePath.find('.') != string::npos;
	}

	static string RemoveExtension(const string& fileName)
	{
		string::size_type const p(fileName.find_last_of('.'));
		return fileName.substr(0, p);
	}

	static string FileNameWithoutExtension(const string& filePath)
	{
		return RemoveExtension(FileName(filePath));
	}

	static string Directory(const string& filePath)
	{
		return filePath.substr(0, filePath.find_last_of("/\\"));
	}

	static string Extension(const string& filePath)
	{
		string::size_type const p(filePath.find_last_of('.'));
		return filePath.substr(p + 1, 3);
	}

	static void CreateDirectoryIfNotExist(const string& dirPath)
	{
		struct stat st;
		if (stat(dirPath.c_str(), &st) != 0)
			mkdir(dirPath.c_str(), 0777);
	}
};