#pragma once
#include "Utils/Types.h"
#include "Utils/Action.h"
using namespace std;

struct ProblemArguments
{
	int Dimension = 2;
	string RHSCode = "sine";
	double Kappa1 = 1;
	double Kappa2 = 1;
	double AnisotropyRatio = 1;
	string Partition = "chiasmus";
};

struct DiscretizationArguments
{
	BigNumber N = 16;
	string Method = "hho";
	string MeshCode = "cart";
	string MeshFilePath = "";
	string Stabilization = "hho";
	string BasisCode = "legendre";
	int PolyDegree = 1;
	bool UsePolynomialSpaceQ = false;
	int PenalizationCoefficient = -1;
	bool StaticCondensation = true;
};

struct MultigridArguments
{
	int Levels = 0;
	int ProlongationCode = 1;
	int MatrixMaxSizeForCoarsestLevel = 1000;
	int WLoops = 1;
	int CellInterpolationCode = 1;
	bool UseGalerkinOperator = false;
	string PreSmootherCode = "bgs";
	string PostSmootherCode = "rbgs";
	int PreSmoothingIterations = 1;
	int PostSmoothingIterations = 1;
	string CoarseningStgyCode = "default";
	CoarseningStrategy CoarseningStgy;
};

struct SolverArguments
{
	string SolverCode = "lu";
	string InitialGuessCode = "0";
	double Tolerance = 1e-8;
	int MaxIterations = 200;
	MultigridArguments MG;
};

struct ProgramArguments
{
	string ActionCodes = "sr";
	Action Actions;
	ProblemArguments Problem;
	DiscretizationArguments Discretization;
	SolverArguments Solver;
	string OutputDirectory = ".";
};