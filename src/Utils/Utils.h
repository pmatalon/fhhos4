#pragma once
#include <functional>
#include "Types.h"
#include "../Mesh/Point.h"
using namespace std;

using RefFunction = function<double(RefPoint)>;
using DomFunction = function<double(DomPoint)>;

class Utils
{
public:
	static int Binomial(int n, int p)
	{
		if (p != 0 && n != p)
			return Binomial(n - 1, p) + Binomial(n - 1, p - 1);
		return 1;
	}

	static int Factorial(int n)
	{
		int f = 1;
		for (int i = 2; i <= n; i++)
			f *= i;
		return f;
	}

	static string MatrixInfo(const SparseMatrix& M, string name)
	{
		double density = (double)M.nonZeros() / (M.rows() * M.cols()) * 100;
		double roundedDensity = (int)(density * 10.0) / 10.0;
		//int roundedDensity = ceil(density);
		return "size(" + name + ")=" + to_string(M.rows()) + "x" + to_string(M.cols()) + ", \tnnz(" + name + ")=" + to_string(M.nonZeros()) + ", \tdensity(" + name + ")=" + to_string(density) + "%";
	}

	static bool FileExists(string filename)
	{
		ifstream ifile(filename);
		return ifile.good();
	}

	static string BeginRed;
	static string BeginGreen;
	static string BeginYellow;
	static string EndColor;

	static void FatalError(string msg)
	{
		cout << Utils::BeginRed << "Error: " << msg << Utils::EndColor << endl;
		cout << "------------------------- FAILURE -------------------------" << endl;
		exit(EXIT_FAILURE);
	}

};
string Utils::BeginRed    = "\033[1;31m";
string Utils::BeginGreen  = "\033[1;32m";
string Utils::BeginYellow = "\033[1;33m";
string Utils::EndColor    = "\033[0m\n";