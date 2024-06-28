#pragma once
#include "Utils/Types.h"
using namespace std;

struct ProblemArguments
{
	EquationType Equation = EquationType::Diffusion;
	int Dimension = -1;
	string GeoCode = "square";
	string TestCaseCode = "";
	string SourceCode = "";
	double HeterogeneityRatio = 1;
	double AnisotropyRatio = 1;
	double AnisotropyAngle = 0; // in radians
	bool ComputeNormalDerivative = false;
	string BCCode = "";
	string Scheme = "g";
};

struct DiscretizationArguments
{
	BigNumber N = 16;
	BigNumber Ny = 0;
	BigNumber Nz = 0;
	string Mesher = "default";
	string MeshCode = "default";
	double Stretch = 0.5;
	string PolyMeshInitialMesh = "cart";
	FaceCoarseningStrategy PolyMeshFaceCoarseningStgy = FaceCoarseningStrategy::InterfaceCollapsing;
	FaceCollapsing PolyMeshBoundaryFaceCollapsing = FaceCollapsing::OnlyCollinear;
	int PolyMeshNAggregPasses = 1;

	string Method = "hho";
	string Stabilization = "hho";
	string BiharStabilization = "hdg";
	string ElemBasisCode = "";
	string FaceBasisCode = "";
	int OrthogonalizeElemBasesCode = -1;
	int OrthogonalizeFaceBasesCode = -1;
	int PolyDegree = 2;
	int RelativeCellPolyDegree = 0;
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
	FaceCollapsing BoundaryFaceCollapsing = FaceCollapsing::OnlyCollinear;
	double FaceCollapsingCoplanarityTolerance = 1e-12;
	ReEntrantCornerMgmt ReEntrantCornerManagement = ReEntrantCornerMgmt::Disabled;
	bool ManageAnisotropy = true; // used only in UAMG
};

struct SolverArguments
{
	string SolverCode = "default";
	string PreconditionerCode = "default";
	string InitialGuessCode = "0";
	StoppingCriteria StoppingCrit = StoppingCriteria::NormalizedResidual;
	double Tolerance = 1e-8;
	double Tolerance2 = 1e-8;
	double StagnationConvRate = 0.90;
	int MaxIterations = 200;
	bool PrintIterationResults = true;
	double RelaxationParameter = 1;
	int BlockSize = -1;
	int Restart = 0;

	string BiHarmonicSolverCode = "fcg";
	string BiHarmonicPreconditionerCode = "s";
	string BiHarmonicPrecSolverCode = "bicgstab";
	double BiHarmonicPrecSolverTol = 0;
	double BiHarmonicPrecSolverMaxIter = 1000;
	bool ComputeIterL2Error = false;
	MultigridArguments MG;
	bool BiHarReconstructBoundary = false;
	int NeighbourhoodDepth = 1;
	int PatchSize = 3;
};

struct ExportArguments
{
	bool LinearSystem = false;
	bool AssemblyTermMatrices = false;
	bool MeshToMatlab = false;
	bool MultigridComponents = false;
	bool MultigridIterationVectors = false;
	bool SolutionVectors = false;
	bool SolutionToGMSH = false;
	bool ExactSolutionToGMSH = false;
	bool ErrorToGMSH = false;
	bool AbsErrorToGMSH = false;
	bool SourceToGMSH = false;
	bool Iterations = false;
	bool IterationResiduals = false;
	bool IterationL2Errors = false;
	string ValueSeparator = ","; // ",\n";
	double VisuTolerance = 1e-3;
	int VisuMaxRefinements = 6;
};

struct ActionsArguments
{
	bool SolveLinearSystem = true;
	bool ComputeErrors = true;
	bool UseCache = true;
	bool LogAssembly = true;
	bool AssembleRightHandSide = true;
	bool InitReferenceShapes = true;
	bool UnitTests = false;
	bool GMSHLogEnabled = false;
	bool PrintDebug = false;
	int Option1 = 0;
	int Option2 = 0;
	double DoubleParam1 = 0.0;
	ExportArguments Export;
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