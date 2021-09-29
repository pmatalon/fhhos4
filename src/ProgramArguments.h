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
	BigNumber Ny = 0;
	BigNumber Nz = 0;
	string Method = "hho";
	string Mesher = "gmsh";
	string MeshCode = "default";
	double Stretch = 0.5;
	string Stabilization = "hho";
	string ElemBasisCode = "";
	string FaceBasisCode = "";
	int OrthogonalizeElemBasesCode = -1;
	int OrthogonalizeFaceBasesCode = -1;
	int PolyDegree = 2;
	bool UsePolynomialSpaceQ = false;
	int PenalizationCoefficient = -1;
	bool StaticCondensation = true;
};

struct MultigridArguments
{
	int Levels = 0;
	int ProlongationCode = 0;
	GMG_H_Prolongation GMG_H_Prolong = GMG_H_Prolongation::Default;
	GMG_P_Prolongation GMG_P_Prolong = GMG_P_Prolongation::Injection;
	GMG_P_Restriction GMG_P_Restrict = GMG_P_Restriction::RemoveHigherOrders;
	int FaceProlongationCode = 0;
	int CoarseningProlongationCode = 0;
	int NSubtriangulationsForApproxL2Proj = 1;
	UAMGFaceProlongation UAMGFaceProlong = UAMGFaceProlongation::BoundaryAggregatesInteriorAverage;
	UAMGProlongation UAMGCoarseningProlong = UAMGProlongation::FaceProlongation;
	UAMGProlongation UAMGMultigridProlong = UAMGProlongation::ReconstructTraceOrInject;
	int NumberOfMeshes = 0;
	int MatrixMaxSizeForCoarsestLevel = 1000;
	char CycleLetter = 'V';
	int WLoops = 1;
	bool UseHigherOrderReconstruction = true;
	bool UseHeterogeneousWeighting = true;
	bool UseGalerkinOperator = false;
	string PreSmootherCode = "bgs";
	string PostSmootherCode = "rbgs";
	int PreSmoothingIterations = 1;
	int PostSmoothingIterations = 1;
	int CoarseLevelChangeSmoothingCoeff = 0;
	char CoarseLevelChangeSmoothingOperator = '+';
	HP_CoarsStgy HP_CS = HP_CoarsStgy::H_only;
	H_CoarsStgy H_CS = H_CoarsStgy::None;
	P_CoarsStgy P_CS = P_CoarsStgy::Minus2;
	FaceCoarseningStrategy FaceCoarseningStgy = FaceCoarseningStrategy::InterfaceCollapsing;
	double CoarseningFactor = 0;
	BigNumber CoarseN = 2;
	string CoarseSolverCode = "ch";
	FaceCollapsing BoundaryFaceCollapsing = FaceCollapsing::Max;
	ReEntrantCornerMgmt ReEntrantCornerManagement = ReEntrantCornerMgmt::Disabled;
	bool ManageAnisotropy = true; // used only in UAMG
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
	bool InitReferenceShapes = true;
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