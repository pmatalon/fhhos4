#pragma once
#include <fstream>
#include <functional>
#include "Types.h"
#include "../Geometry/Point.h"
using namespace std;

#ifndef ROOT_PATH
#define ROOT_PATH "./"
#endif // !ROOT_PATH

using RefFunction = function<double(RefPoint)>;
using DomFunction = function<double(DomPoint)>;

class Utils
{
public:
	static string RootPath()
	{
		return ROOT_PATH;
	}

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

	static bool IsPowerOf2(BigNumber n)
	{
		return (n & (n - 1)) == 0;
	}

	static string MatrixInfo(const SparseMatrix& M, string name)
	{
		double density = (double)M.nonZeros() / (M.rows() * M.cols()) * 100;
		//double roundedDensity = (int)(density * 10.0) / 10.0;
		//int roundedDensity = ceil(density);
		return "size(" + name + ")=" + to_string(M.rows()) + "x" + to_string(M.cols()) + ", \tnnz(" + name + ")=" + to_string(M.nonZeros()) + ", \tdensity(" + name + ")=" + to_string(density) + "%";
	}

	static void Empty(DenseMatrix& M)
	{
		M = DenseMatrix(0, 0);
	}
	static void Empty(SparseMatrix& M)
	{
		M = SparseMatrix(0, 0);
	}
	static void Empty(Vector& v)
	{
		v = Vector(0);
	}

	static size_t MemoryUsage(const SparseMatrix& M)
	{
		return M.nonZeros() * (sizeof(double) + 2* sizeof(Eigen::Index));
	}
	static size_t MemoryUsage(const Vector& v)
	{
		return v.rows() * sizeof(double);
	}

	static string MemoryString(size_t bytes)
	{
		stringstream ss;
		ss.precision(2);

		double kilos = bytes / 1024.0;
		double megs = kilos / 1024.0;
		double gigs = megs / 1024.0;
		if (gigs >= 1)
			ss << gigs << "GB";
		else if (megs >= 1)
			ss << megs << "MB";
		else
			ss << kilos << "KB";
		return ss.str();
	}

	static bool FileExists(string filename)
	{
		ifstream ifile(filename);
		return ifile.good();
	}

	template <typename T>
	static inline vector<T> Join(const vector<T>& A, const vector<T>& B)
	{
		vector<T> AB;
		AB.reserve(A.size() + B.size()); // preallocate memory
		AB.insert(AB.end(), A.begin(), A.end());
		AB.insert(AB.end(), B.begin(), B.end());
		return AB;
	}

	template <typename T>
	static inline vector<T> SymmetricDifference(const vector<T>& A, const vector<T>& B)
	{
		vector<T> v1 = A;
		vector<T> v2 = B;

		sort(v1.begin(), v1.end());
		sort(v2.begin(), v2.end());

		vector<T> symDifference;

		set_symmetric_difference(v1.begin(), v1.end(),
								 v2.begin(), v2.end(),
								 std::back_inserter(symDifference));
		return symDifference;
	}

	template <typename T>
	static inline bool HasDuplicates(vector<T> v)
	{
		sort(v.begin(), v.end());
		auto it = unique(v.begin(), v.end());
		bool wasUnique = (it == v.end());
		return !wasUnique;
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

	static vector<string> Explode(const string& stringList, char separator)
	{
		stringstream ss(stringList);
		vector<string> result;
		while (ss.good())
		{
			string substr;
			getline(ss, substr, separator);
			result.push_back(substr);
		}
		return result;
	}

	static bool IsPredefinedGeometry(const string& geo)
	{
		return geo.compare("segment") == 0
			|| geo.compare("square") == 0
			|| geo.compare("square4quadrants") == 0
			|| geo.compare("cube") == 0;
	}

	static bool IsRefinementStrategy(CoarseningStrategy stgy)
	{
		return stgy == CoarseningStrategy::BeyRefinement || stgy == CoarseningStrategy::GMSHSplittingRefinement;
	}

	static bool BuildsNestedMeshHierarchy(CoarseningStrategy stgy)
	{
		return IsRefinementStrategy(stgy) || stgy == CoarseningStrategy::StandardCoarsening;
	}

	static bool RequiresNestedHierarchy(Prolongation p)
	{
		return p != Prolongation::CellInterp_L2proj_Trace;
	}

	static double Eps;
	static double NumericalZero;

	static string BeginRed;
	static string BeginGreen;
	static string BeginYellow;
	static string EndColor;

	static void Error(string msg)
	{
		cout << Utils::BeginRed << "Error: " << msg << Utils::EndColor << endl;
	}
	static void FatalError(string msg)
	{
		cout << Utils::BeginRed << "Error: " << msg << Utils::EndColor << endl;
		cout << "------------------------- FAILURE -------------------------" << endl;
		exit(EXIT_FAILURE);
	}
	static void Warning(string msg)
	{
		Warning(cout, msg);
	}
	static void Warning(ostream& os, string msg)
	{
		os << Utils::BeginYellow << "Warning: " << msg << Utils::EndColor << endl;
	}

};

double Utils::Eps = 1e-4;
double Utils::NumericalZero = 1e-12;

string Utils::BeginRed    = "\033[1;31m";
string Utils::BeginGreen  = "\033[1;32m";
string Utils::BeginYellow = "\033[1;33m";
string Utils::EndColor    = "\033[0m";