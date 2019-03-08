#pragma once
#include <vector>
#include <Eigen/Sparse>

class NonZeroCoefficients
{
private:
	vector<Eigen::Triplet<double>> coefficients;
	Eigen::MatrixXf Done;
public:
	NonZeroCoefficients(BigNumber nnzApproximate)
	{
		this->coefficients.reserve(nnzApproximate);
	}
	NonZeroCoefficients(BigNumber nRows, BigNumber nCols, BigNumber nnzApproximate) : Done(nRows, nCols)
	{
		this->coefficients.reserve(nnzApproximate);
	}
	void Add(BigNumber i, BigNumber j, double value)
	{
		if (Done.rows() > 0)
		{
			if (Done(i, j) == 1)
				return;
			Done(i, j) = 1;
		}

		if (value != 0)
			this->coefficients.push_back(Eigen::Triplet<double>(i, j, value));
	}
	void Fill(Eigen::SparseMatrix<double> &m)
	{
		m.setFromTriplets(this->coefficients.begin(), this->coefficients.end());
	}
};