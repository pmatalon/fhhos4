#pragma once
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <vector>
#include <unistd.h>
using namespace std;

// Directory of the data (meshes/) relative to the directory of the executable: share/fhhos4 in the installation, and in
// the build tree, where CMake links its meshes/ to data/meshes.
#ifndef DATA_DIR_FROM_BIN
#define DATA_DIR_FROM_BIN "../share/fhhos4"
#endif // !DATA_DIR_FROM_BIN

class FileSystem
{
public:
	// Directory of the data (meshes/): $FHHOS4_DATA_DIR if set, DATA_DIR_FROM_BIN from the directory of the executable
	// otherwise.
	static string DataDirectory()
	{
		const char* dir = getenv("FHHOS4_DATA_DIR");
		if (dir && *dir)
			return dir;
		error_code ec;
		filesystem::path exe = filesystem::read_symlink("/proc/self/exe", ec); // Linux
		if (ec)
			return DATA_DIR_FROM_BIN; // from the working directory
		return (exe.parent_path() / DATA_DIR_FROM_BIN).lexically_normal().string();
	}

	// Message of the error "data file not found", after a search at the given paths, with how to give the data directory.
	// Every failed search of a data file (e.g. a mesh) ends with it.
	static string DataFileNotFound(const string& file, const vector<string>& searchedPaths)
	{
		string msg = "File not found: " + file + ". Searched:";
		for (const string& path : searchedPaths)
			msg += "\n    " + path;
		const char* dir = getenv("FHHOS4_DATA_DIR");
		if (dir && *dir)
			msg += "\nFHHOS4_DATA_DIR is set to " + string(dir) + ": it must be the directory data of the sources of fhhos4:";
		else
			msg += "\nThe data of fhhos4 (data/ in its sources) is searched in " + DataDirectory() + ", next to the directory "
				"of the executable: in its build directory or its installation. Elsewhere, give the directory data of the sources:";
		msg += "\n    export FHHOS4_DATA_DIR=<path to fhhos4>/data";
		return msg;
	}

	// Cache (GMSH meshes): $XDG_CACHE_HOME/fhhos4, by default ~/.cache/fhhos4.
	static string CacheDirectory()
	{
		const char* cache = getenv("XDG_CACHE_HOME");
		if (cache && *cache)
			return string(cache) + "/fhhos4";
		const char* home = getenv("HOME");
		if (home && *home)
			return string(home) + "/.cache/fhhos4";
		return (filesystem::temp_directory_path() / "fhhos4_cache").string();
	}

	// Path of a temporary file, unique to the process: runs at the same time don't share it. Removed at the end of the
	// process at the latest.
	static string TemporaryFile(const string& name)
	{
		string path = (filesystem::temp_directory_path() / ("fhhos4_" + to_string(getpid()) + "_" + name)).string();
		TemporaryFiles().push_back(path);
		return path;
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

	// With its parent directories
	static void CreateDirectoryIfNotExist(const string& dirPath)
	{
		error_code ec;
		filesystem::create_directories(dirPath, ec);
	}

private:
	struct TemporaryFileList : vector<string>
	{
		~TemporaryFileList()
		{
			for (const string& path : *this)
				remove(path.c_str());
		}
	};
	// Destroyed at the end of the process (static)
	static TemporaryFileList& TemporaryFiles()
	{
		static TemporaryFileList files;
		return files;
	}
};