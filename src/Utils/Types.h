#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>

typedef long unsigned int BigNumber;

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

using DenseMatrix = Eigen::MatrixXd;

using Vector = Eigen::VectorXd;

template<int Dim>
using DimVector = Eigen::Matrix<double, Dim, 1>;

template<int Dim>
using DimMatrix = Eigen::Matrix<double, Dim, Dim>;