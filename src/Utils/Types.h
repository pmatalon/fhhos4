#pragma once
#include <iomanip>
#include <cstdio>
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <Eigen/Sparse>
#include <Eigen/Dense>

// BigNumber is signed in both cases: the code encodes states as negative numbers (e.g. -1: not aggregated yet)
#ifdef SMALL_INDEX
using BigNumber = int;
using SparseMatrixIndex = int;
#else
using BigNumber = intmax_t;
using SparseMatrixIndex = intmax_t; // in Eigen, default value is int, but it is not sufficient for big matrices
#endif // SMALL_INDEX

using RowMajorSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, SparseMatrixIndex>;
using ColMajorSparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, SparseMatrixIndex>;

// Row-major because the smoothers (Gauss-Seidel, SOR, block Jacobi) and the AMG setups iterate over rows.
// Eigen's sparse direct solvers only support column-major matrices: they must be given a column-major
// copy (see EigenSparseLU, EigenSparseCholesky and IterativeSolver).
using SparseMatrix = RowMajorSparseMatrix;


using DenseMatrix = Eigen::MatrixXd;

using Vector = Eigen::VectorXd;

template<int Dim>
using DimVector = Eigen::Matrix<double, Dim, 1>;

template<int Dim>
using DimMatrix = Eigen::Matrix<double, Dim, Dim>;

using Flops = std::size_t;
using MFlops = double;

enum class EquationType : unsigned
{
	Diffusion,
	BiHarmonic,
	BiHarmonicDD
};

enum class PbBoundaryConditions : unsigned
{
	FullDirichlet,
	MixedDirichletNeumann,
	FullNeumann
};

enum class BoundaryConditionType : unsigned
{
	NotOnBoundary = 0,
	Dirichlet = 1,
	Neumann = 2
};

enum class StoppingCriteria : unsigned
{
	NormalizedResidual,
	Stagnation,
	MaxIterations,
};

// Flexible Krylov method: FCG for symmetric positive definite matrices, GCR for the others (e.g. in the K-cycle)
enum class KrylovMethod : unsigned
{
	FCG,
	GCR
};

// Multigrid cycle (-cycle)
enum class CycleType : unsigned
{
	V,
	W, // Multigrid::WLoops loops
	K  // each coarse level is solved by 2 iterations of a flexible Krylov method (Multigrid::KCycleKrylov)
};

// How the numbers of smoothing iterations change from a level to the next coarser one (-cycle V,1,1,+1)
enum class SmoothingChangeOperator : unsigned
{
	Add,
	Subtract,
	Multiply
};

// The smoothers of the multigrid methods (-smoothers; their codes: SmootherFactory). The block versions have one block
// per face (or cell); "relaxed": with the relaxation parameter (-relax, default 1).
enum class SmootherType : unsigned
{
	Default, // post-smoother: the one that makes the cycle symmetric with the pre-smoother (-smoothers with one code)
	Jacobi, // relaxed
	DampedJacobi, // damping factor 2/3
	BlockJacobi, // relaxed
	DampedBlockJacobi, // damping factor 2/3
	GaussSeidel,
	BackwardGaussSeidel,
	SymmetricGaussSeidel, // forward then backward
	AlternatingGaussSeidel, // forward and backward alternately, starting with forward
	BackwardAlternatingGaussSeidel, // the same, starting with backward
	SOR, // relaxed Gauss-Seidel
	BackwardSOR,
	BlockGaussSeidel, // relaxed (block SOR)
	BackwardBlockGaussSeidel,
	SymmetricBlockGaussSeidel,
	AlternatingBlockGaussSeidel,
	BackwardAlternatingBlockGaussSeidel,
	HybridBlockGaussSeidel, // Gauss-Seidel in the rows of each thread, Jacobi between the threads
	HybridBackwardBlockGaussSeidel,
	ChebyshevFirstKind, // preconditioned by Jacobi
	BlockChebyshevFirstKind, // preconditioned by block Jacobi
	ChebyshevFourthKind,
	BlockChebyshevFourthKind
};

enum class CoarseningType : unsigned
{
	H,
	P,
	HP
};

enum class HP_CoarsStgy : unsigned
{
	H_only,
	P_only,
	P_then_H,
	H_then_P, // the number of meshes must be fixed
	HP_then_H,
	HP_then_P, // the number of meshes must be fixed
	P_then_HP, // the number of meshes must be fixed
	Alternate
};

enum class H_CoarsStgy : unsigned
{
	None,
	StandardCoarsening,
	AgglomerationCoarseningByVertexRemoval,
	AgglomerationCoarseningByClosestCenter,
	AgglomerationCoarseningByClosestFace,
	AgglomerationCoarseningByLargestInterface,
	AgglomerationCoarseningBySeedPoints,
	AgglomerationCoarseningByMostCoplanarFaces,
	AgglomerationCoarseningByFaceNeighbours,
	AgglomerationCoarseningByVertexNeighbours,
	MultipleAgglomerationCoarseningByFaceNeighbours,
	IndependentRemeshing,
	FaceCoarsening,
	GMSHSplittingRefinement,
	BeyRefinement,
	DoublePairwiseAggregation,
	MultiplePairwiseAggregation
};

enum class P_CoarsStgy : unsigned
{
	Minus1,
	Minus2,
	DivideBy2,
	DirectToLow
};

enum class FaceCoarseningStrategy : unsigned
{
	None,
	InterfaceCollapsing,
	InterfaceCollapsingAndTryAggregInteriorToInterfaces
};

enum class GMG_H_Prolongation : unsigned
{
	Default = 0,
	CellInterp_Trace = 1,
	CellInterp_Inject_Trace = 6,
	CellInterp_ExactL2proj_Trace = 7,
	CellInterp_ApproxL2proj_Trace = 8,
	CellInterp_FinerApproxL2proj_Trace = 9,
	CellInterp_InjectAndTrace = 2,
	CellInterp_Inject_Adjoint = 3,
	Wildey = 4,
	FaceInject = 5
};

enum class GMG_P_Prolongation : unsigned
{
	Injection = 1,
	H_Prolongation = 2
};

enum class GMG_P_Restriction : unsigned
{
	RemoveHigherOrders = 1,
	P_Transpose = 2
};

// The prolongations of U-AMG: the coarsening prolongation of each pass (-coarsening-prolong) and the multigrid
// prolongation of a level (-prolong). Their use: development (experiments), the 2023 AMG paper
// (reproducibility/2023_AMG_for_hybrid_methods.md), the library fhhos4 (fhhos4::Prolongation).
enum class UAMGProlongation : unsigned
{
	ReconstructionTrace = 1,                  // development use
	ChainedCoarseningProlongations = 2,       // used in the paper: the multigrid prolongation (-prolong 2, the default),
	                                          // the coarsening prolongations chained
	FaceProlongation = 3,                     // used in the paper (Q_F)
	FaceProlongationAndInteriorSmoothing = 4, // used in the paper (Q_F^smooth)
	ReconstructTraceOrInject = 5,             // used in the paper (P_F^(0))
	ReconstructSmoothedTraceOrInject = 6,     // used in the paper (P_F, the program's default); library
	                                          // (Prolongation::Reconstruction)
	FindInteriorThatReconstructs = 7,         // development use
	HighOrder = 8,                            // development use
	ReconstructionTranspose2Steps = 9,        // development use
	// Coarsening prolongation only: injection on the kept faces, local solve of their rows of the operator on the
	// removed faces (FaceProlongationAndInteriorSmoothing with omega = 1). See CoarseningPass::ComputeProlongation().
	HarmonicExtensionOrInject = 10            // library (Prolongation::HarmonicExtension, its default)
};

enum class UAMGFaceProlongation : unsigned
{
	BoundaryAggregatesInteriorAverage = 1,
	BoundaryAggregatesInteriorZero = 2,
	FaceAggregates = 3,
};

enum class FaceCollapsingStatus : unsigned
{
	Ok,
	NotEnoughFaces,
	GeometricErosion,
	InterfaceHasHoles,
	ElementFullDegeneration,
	ElementPartialDegeneration,
	CrossedPolygon,
	OneElementEmbeddedInConvexHullOfTheOther
};

enum class FaceCollapsing : unsigned
{
	Disabled,
	OnlyCollinear,
	ByPairs,
	Max
};

enum class ReEntrantCornerMgmt : unsigned
{
	Disabled,
	AgglomerateFirst
};

enum PolygonalTriangulation
{
	Barycentric,
	OneVertex
};



enum class ToleranceStrategy : unsigned
{
	Fixed,
	DynamicFixedStep,
	DynamicVariableStep
};