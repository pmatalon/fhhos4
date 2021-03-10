#pragma once
#include "Utils/Types.h"
using namespace std;

struct ProblemArguments
{
	int Dimension = -1;
	string GeoCode = "square";
	string TestCaseCode = "";
	string SourceCode = "";
	double HeterogeneityRatio = 1;
	double AnisotropyRatio = 1;
	double AnisotropyAngle = 0; // in radians
	string BCCode = "d";
};

struct DiscretizationArguments
{
	BigNumber N = 16;
	BigNumber Ny = -1;
	BigNumber Nz = -1;
	string Method = "hho";
	string Mesher = "gmsh";
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
	int ProlongationCode = 0;
	Prolongation GMGProlong = Prolongation::Default;
	int FaceProlongationCode = 0;
	CAMGFaceProlongation CAMGFaceProlong = CAMGFaceProlongation::BoundaryAggregatesInteriorAverage;
	CAMGProlongation CAMGProlong = CAMGProlongation::ReconstructTraceOrInject;
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
	CoarseningStrategy CoarseningStgy = CoarseningStrategy::None;
	FaceCoarseningStrategy FaceCoarseningStgy = FaceCoarseningStrategy::InterfaceCollapsing;
	double CoarseningFactor = 0;
	BigNumber CoarseN = 2;
	string CoarseSolverCode = "ch";
	FaceCollapsing BoundaryFaceCollapsing = FaceCollapsing::Max;
	ReEntrantCornerMgmt ReEntrantCornerManagement = ReEntrantCornerMgmt::Disabled;
	bool ManageAnisotropy = true; // used only in CAMG
};

struct SolverArguments
{
	string SolverCode = "default";
	string InitialGuessCode = "0";
	double Tolerance = 1e-8;
	int MaxIterations = 200;
	bool PrintIterationResults = true;
	double RelaxationParameter = 1;
	int BlockSize = -1;
	MultigridArguments MG;
};

struct ActionsArguments
{
	bool SolveLinearSystem = true;
	bool UseCache = true;
	bool ExportLinearSystem = false;
	bool ExportAssemblyTermMatrices = false;
	bool ExportMeshToMatlab = false;
	bool ExportMultigridComponents = false;
	bool ExportMultigridIterationVectors = false;
	bool ExportSolutionVectors = false;
	bool ExportSolutionToGMSH = false;
	bool ExportErrorToGMSH = false;
	bool ExportSourceToGMSH = false;
	bool LogAssembly = true;
	bool AssembleRightHandSide = true;
	bool UnitTests = false;
	bool GMSHLogEnabled = false;
};

struct ProgramArguments
{
	string ActionCodes = "sr";
	ActionsArguments Actions;
	ProblemArguments Problem;
	DiscretizationArguments Discretization;
	SolverArguments Solver;
	string OutputDirectory = "./out";
};