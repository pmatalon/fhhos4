#pragma once
#include <fstream>
#include <functional>
#include "Types.h"
#include "ParallelLoop.h"
#include "../Geometry/Point.h"
#include "../ProgramArguments.h"
using namespace std;

using RefFunction = function<double(const RefPoint&)>;
using DomFunction = function<double(const DomPoint&)>;

class Utils
{
public:

	static ProgramArguments ProgramArgs;

	//-------------------//
	// Util DomFunctions //
	//-------------------//

	static double ConstantFunctionZero(const DomPoint& p)
	{
		return 0;
	}
	static double ConstantFunctionOne(const DomPoint& p)
	{
		return 1;
	}

	//-------------//
	// Computation //
	//-------------//

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

	static bool IsInDisk(const DomPoint& c, double r, const DomPoint& p)
	{
		if (pow(p.X - c.X, 2) + pow(p.Y - c.Y, 2) <= r)
			return true;
		return false;
	}

	//--------------//
	//   Matrices   //
	//--------------//

	static string MatrixInfo(const SparseMatrix& M, string name)
	{
		double density = (double)M.nonZeros() / (M.rows() * M.cols()) * 100;
		//double roundedDensity = (int)(density * 10.0) / 10.0;
		//int roundedDensity = ceil(density);
		return "size(" + name + ")=" + to_string(M.rows()) + "x" + to_string(M.cols()) + ", \tnnz(" + name + ")=" + to_string(M.nonZeros()) + ", \tdensity(" + name + ")=" + to_string(density) + "%";
	}

	static SparseMatrix InvertBlockDiagMatrix(const SparseMatrix& M, int blockSize)
	{
		BigNumber nElements = M.rows() / blockSize;
		NumberParallelLoop<> parallelLoop(nElements);
		parallelLoop.ReserveChunkCoeffsSize(blockSize * blockSize);
		parallelLoop.Execute([&M, blockSize](BigNumber i, ParallelChunk<CoeffsChunk>* chunk)
			{
				DenseMatrix block = M.block(i * blockSize, i * blockSize, blockSize, blockSize);
				chunk->Results.Coeffs.Add(i * blockSize, i * blockSize, block.inverse());
			});
		SparseMatrix invM(M.rows(), M.cols());
		parallelLoop.Fill(invM);
		return invM;
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

	static double Cond(const DenseMatrix& A)
	{
		Eigen::EigenSolver<DenseMatrix> es(A);
		Eigen::VectorXcd eigenvalues = es.eigenvalues();
		auto realEigenvalues = eigenvalues.real();
		return realEigenvalues.maxCoeff() / realEigenvalues.minCoeff();
	}

	//------------//
	//   Memory   //
	//------------//

	static size_t MemoryUsage(const SparseMatrix& M)
	{
		return M.nonZeros() * (sizeof(double) + 2* sizeof(SparseMatrix::StorageIndex));
	}
	static size_t MemoryUsage(const Vector& v)
	{
		return v.rows() * sizeof(double);
	}
	static size_t SparseMatrixMemoryUsage(size_t nonZeroes)
	{
		return nonZeroes * (sizeof(double) + 2 * sizeof(SparseMatrix::StorageIndex));
	}
	static size_t DenseMatrixMemoryUsage(BigNumber rows, BigNumber cols)
	{
		return cols * rows * sizeof(double) + sizeof(DenseMatrix);
	}
	static size_t VectorMemoryUsage(BigNumber rows)
	{
		return rows * sizeof(double);
	}

	static string MemoryString(size_t bytes)
	{
		stringstream ss;
		ss << fixed;
		ss.precision(0);

		double kilos = bytes / 1024.0;
		double megs = kilos / 1024.0;
		double gigs = megs / 1024.0;
		if (gigs >= 1)
		{
			ss.precision(1);
			ss << gigs << "GB";
		}
		else if (megs >= 1)
			ss << megs << "MB";
		else
			ss << kilos << "KB";
		return ss.str();
	}

	//-------------//
	//   Vectors   //
	//-------------//

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

	//---------------------------//
	//   Coarsening/Refinement   //
	//---------------------------//

	static bool IsPredefinedGeometry(const string& geo)
	{
		return geo.compare("segment") == 0
			|| geo.compare("square") == 0
			|| geo.compare("square4quadrants") == 0
			|| geo.compare("cube") == 0;
	}

	static bool IsRefinementStrategy(H_CoarsStgy stgy)
	{
		return stgy == H_CoarsStgy::BeyRefinement || stgy == H_CoarsStgy::GMSHSplittingRefinement;
	}

	static bool BuildsNestedMeshHierarchy(H_CoarsStgy stgy)
	{
		return IsRefinementStrategy(stgy) || stgy == H_CoarsStgy::StandardCoarsening;
	}

	static bool IsAlgebraic(string solverCode)
	{
		return EndsWith(solverCode, "amg");
	}

	static bool RequiresNestedHierarchy(GMG_H_Prolongation p)
	{
		return p != GMG_H_Prolongation::CellInterp_ExactL2proj_Trace && p != GMG_H_Prolongation::CellInterp_ApproxL2proj_Trace && p != GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace;
	}

	static double Eps;
	static double NumericalZero;

	//-------------//
	//   Logging   //
	//-------------//

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
		assert(false);
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


	static bool Contains(const std::string& str1, const std::string& str2)
	{
		return str1.find(str2) != std::string::npos;
	}
	static bool StartsWith(const std::string& str, const std::string& prefix)
	{
		return str.rfind(prefix, 0) == 0;
	}
	static bool EndsWith(const std::string& str, const std::string& suffix)
	{
		return str.size() >= suffix.size() && 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
	}

};

ProgramArguments Utils::ProgramArgs;

double Utils::Eps = 1e-4;
double Utils::NumericalZero = 1e-12;

string Utils::BeginRed    = "\033[1;31m";
string Utils::BeginGreen  = "\033[1;32m";
string Utils::BeginYellow = "\033[1;33m";
string Utils::EndColor    = "\033[0m";