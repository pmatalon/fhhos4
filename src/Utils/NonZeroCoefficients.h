#pragma once
#include <vector>
#include <Eigen/Sparse>
#include "Utils.h"

class NonZeroCoefficients
{
private:
	vector<Eigen::Triplet<double>> coefficients;
public:
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
	inline void Add(BigNumber i, BigNumber j, double value)
	{
		if (abs(value) > 1e-16)
			this->coefficients.push_back(Eigen::Triplet<double>(i, j, value));
	}

	inline void Add(NonZeroCoefficients &chunk)
	{
		this->coefficients.insert(this->coefficients.end(), chunk.coefficients.begin(), chunk.coefficients.end());
	}

	inline void Add(BigNumber iStart, BigNumber jStart, const DenseMatrix &m)
	{
		for (int i = 0; i < m.rows(); ++i)
		{
			for (int j = 0; j < m.cols(); ++j)
				Add(iStart + i, jStart + j, m(i, j));
		}
	}

	inline void CopyRows(BigNumber iStart, int nRows, const Eigen::SparseMatrix<double, Eigen::RowMajor> &m)
	{
		for (int i = 0; i < nRows; ++i)
		{
			for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(m, iStart + i); it; ++it)
			{
				Add(it.row(), it.col(), it.value());
			}
		}
	}

	inline void Fill(SparseMatrix &m)
	{
		m.setFromTriplets(this->coefficients.begin(), this->coefficients.end());
	}

	friend ostream& operator<<(ostream& os, const NonZeroCoefficients& nnzCoeffs)
	{
		for (Eigen::Triplet<double> t : nnzCoeffs.coefficients)
			os << t.row() << "\t" << t.col() << "\t" << t.value() << endl;
		return os;
	}
};
