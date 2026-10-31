#pragma once
#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace fhhos4_tests
{
	// M ⊗ C: the coefficient (i, j) of M times C, at the rows i C.rows() + r and the columns j C.cols() + s. With C = I_n,
	// the vector-valued version of a scalar block, its components interleaved (fhhos4::Solver::VectorComponents).
	template <typename Matrix>
	Matrix Kron(const Matrix& M, const Eigen::MatrixXd& C)
	{
		std::vector<Eigen::Triplet<double, typename Matrix::StorageIndex>> triplets;
		for (int i = 0; i < M.outerSize(); i++)
			for (typename Matrix::InnerIterator it(M, i); it; ++it)
				for (int r = 0; r < C.rows(); r++)
					for (int s = 0; s < C.cols(); s++)
						if (C(r, s) != 0)
							triplets.emplace_back(it.row() * C.rows() + r, it.col() * C.cols() + s, it.value() * C(r, s));
		Matrix K(M.rows() * C.rows(), M.cols() * C.cols());
		K.setFromTriplets(triplets.begin(), triplets.end());
		return K;
	}

	// The vector-valued version of a scalar block: M ⊗ I_n
	template <int Options, typename Index>
	Eigen::SparseMatrix<double, Options, Index> OnEachComponent(const Eigen::SparseMatrix<double, Options, Index>& M, int n)
	{
		return Kron(M, Eigen::MatrixXd::Identity(n, n));
	}

	// v ⊗ 1_n: the same value on the n components
	inline Eigen::VectorXd OnEachComponent(const Eigen::VectorXd& v, int n)
	{
		return v.replicate(1, n).transpose().reshaped();
	}
}
