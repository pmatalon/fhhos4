#pragma once
#include <Eigen/Sparse>

typedef long unsigned int BigNumber;

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

template<int Dim>
using DimVector = Eigen::Matrix<double, Dim, 1>;

template<int Dim>
using DimMatrix = Eigen::Matrix<double, Dim, Dim>;