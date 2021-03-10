#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>

#ifdef SMALL_INDEX
using BigNumber = int;
using SparseMatrixIndex = int;
#else
using BigNumber = std::size_t;
using SparseMatrixIndex = intmax_t; // in Eigen, default value is int, but it is not sufficient for big matrices
#endif // SMALL_INDEX

using RowMajorSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, SparseMatrixIndex>;
using ColMajorSparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, SparseMatrixIndex>;

using SparseMatrix = RowMajorSparseMatrix;


using DenseMatrix = Eigen::MatrixXd;

using Vector = Eigen::VectorXd;

template<int Dim>
using DimVector = Eigen::Matrix<double, Dim, 1>;

template<int Dim>
using DimMatrix = Eigen::Matrix<double, Dim, Dim>;

enum class BoundaryConditionType : unsigned
{
	NotOnBoundary = 0,
	Dirichlet = 1,
	Neumann = 2
};

enum class CoarseningStrategy : unsigned
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
	IndependentRemeshing,
	FaceCoarsening,
	GMSHSplittingRefinement,
	BeyRefinement,
	DoublePairwiseAggregation,
	MultiplePairwiseAggregation
};

enum class FaceCoarseningStrategy : unsigned
{
	None,
	InterfaceCollapsing,
	InterfaceCollapsingAndTryAggregInteriorToInterfaces
};

enum class Prolongation : unsigned
{
	Default = 0,
	CellInterp_Trace = 1,
	CellInterp_Inject_Trace = 6,
	CellInterp_L2proj_Trace = 7,
	CellInterp_ApproxL2proj_Trace = 8,
	CellInterp_FinerApproxL2proj_Trace = 9,
	CellInterp_InjectAndTrace = 2,
	CellInterp_Inject_Adjoint = 3,
	Wildey = 4,
	FaceInject = 5
};

enum class CAMGProlongation : unsigned
{
	ReconstructionTrace1Step = 1,
	ReconstructionTrace2Steps = 2,
	FaceProlongation = 3,
	ReconstructionTranspose2Steps = 4,
	ReconstructTraceOrInject = 5,
	ReconstructSmoothedTraceOrInject = 6,
	FindInteriorThatReconstructs = 7,
	HighOrder = 8
};

enum class CAMGFaceProlongation : unsigned
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