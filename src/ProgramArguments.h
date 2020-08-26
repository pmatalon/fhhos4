#pragma once
#include "Utils/Types.h"
#include "Utils/Action.h"
using namespace std;

struct ProblemArguments
{
	int Dimension = -1;
	string GeoCode = "square";
	string TestCaseCode = "sine";
	double HeterogeneityRatio = 1;
	double AnisotropyRatio = 1;
	double AnisotropyAngleInDegrees = 0;
	string Partition = "chiasmus";
	string BCCode = "d";
};

struct DiscretizationArguments
{
	BigNumber N = 16;
	BigNumber Ny = -1;
	BigNumber Nz = -1;
	string Method = "hho";
	string Mesher = "inhouse";
	string MeshCode = "default";
	double Stretch = 0.5;
	string Stabilization = "hho";
	string BasisCode = "legendre";
	int PolyDegree = 2;
	bool UsePolynomialSpaceQ = false;
	int PenalizationCoefficient = -1;
	bool StaticCondensation = true;
};

struct MultigridArguments
{
	int Levels = 0;
	Prolongation ProlongationCode = Prolongation::CellInterp_L2proj_Trace;
	int MatrixMaxSizeForCoarsestLevel = 1000;
	char CycleLetter = 'V';
	int WLoops = 1;
	int CellInterpolationCode = 1;
	string WeightCode = "k";
	bool UseGalerkinOperator = false;
	string PreSmootherCode = "bgs";
	string PostSmootherCode = "rbgs";
	int PreSmoothingIterations = 1;
	int PostSmoothingIterations = 1;
	int CoarseLevelChangeSmoothingCoeff = 0;
	char CoarseLevelChangeSmoothingOperator = '+';
	string CoarseningStgyCode = "default";
	CoarseningStrategy CoarseningStgy;
};

struct SolverArguments
{
	string SolverCode = "default";
	string InitialGuessCode = "0";
	double Tolerance = 1e-8;
	int MaxIterations = 200;
	double RelaxationParameter = 1;
	MultigridArguments MG;
};

struct ProgramArguments
{
	string ActionCodes = "sr";
	Action Actions;
	ProblemArguments Problem;
	DiscretizationArguments Discretization;
	SolverArguments Solver;
	string OutputDirectory = "./out";
};