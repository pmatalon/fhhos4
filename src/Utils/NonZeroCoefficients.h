#pragma once
#include <vector>
#include <Eigen/Sparse>

class NonZeroCoefficients
{
private:
	vector<Eigen::Triplet<double>> coefficients;
	//Eigen::MatrixXf Done;
public:
	int test = 0;
	NonZeroCoefficients(BigNumber nnzApproximate)
	{
		this->coefficients.reserve(nnzApproximate);
	}

	NonZeroCoefficients() {}

	//Reserve()
	/*NonZeroCoefficients(BigNumber nRows, BigNumber nCols, BigNumber nnzApproximate) : Done(nRows, nCols)
	{
		this->coefficients.reserve(nnzApproximate);
	}*/
	void Add(BigNumber i, BigNumber j, double value)
	{
		/*if (Done.rows() > 0)
		{
			if (Done(i, j) == 1)
				return;
			Done(i, j) = 1;
		}*/

		if (abs(value) > 1e-16)
			this->coefficients.push_back(Eigen::Triplet<double>(i, j, value));
	}

	void Add(NonZeroCoefficients &chunk)
	{
		this->coefficients.insert(this->coefficients.end(), chunk.coefficients.begin(), chunk.coefficients.end());
	}

	void Add(BigNumber iStart, BigNumber jStart, const Eigen::MatrixXd &m)
	{
		for (int i = 0; i < m.rows(); ++i)
		{
			for (int j = 0; j < m.cols(); ++j)
				Add(iStart + i, jStart + j, m(i, j));
		}
	}

	void Fill(Eigen::SparseMatrix<double> &m)
	{
		m.setFromTriplets(this->coefficients.begin(), this->coefficients.end());
	}
};
