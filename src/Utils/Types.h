#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>

typedef std::size_t BigNumber;

using SparseMatrixIndex = intmax_t; // in Eigen, default value is int, but it is not sufficient for big matrices

using RowMajorSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, SparseMatrixIndex>;

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
	AgglomerationCoarsening,
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
	BeyRefinement
};

enum class Prolongation : unsigned
{
	Default = 0,
	CellInterp_Trace = 1,
	CellInterp_Inject_Trace = 6,
	CellInterp_L2proj_Trace = 7,
	CellInterp_ApproxL2proj_Trace = 8,
	CellInterp_InjectAndTrace = 2,
	CellInterp_Inject_Adjoint = 3,
	Wildey = 4,
	FaceInject = 5
};

enum class FaceCollapsingStatus : unsigned
{
	Ok,
	NotEnoughFaces,
	InterfaceHasHoles,
	ElementFullDegeneration,
	ElementPartialDegeneration,
	CrossedPolygon,
	OneElementEmbeddedInConvexHullOfTheOther
};