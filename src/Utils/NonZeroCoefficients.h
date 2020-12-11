#pragma once
#include <vector>
#include "Utils.h"

class NonZeroCoefficients
{
protected:
	vector<Eigen::Triplet<double, SparseMatrixIndex>> coefficients;
public:
	NonZeroCoefficients(BigNumber nnzApproximate)
	{
		this->coefficients.reserve(nnzApproximate);
	}

	NonZeroCoefficients() {}

	void Reserve(BigNumber nnzApproximate)
	{
		this->coefficients.reserve(nnzApproximate);
	}

	size_t Size()
	{
		return coefficients.size();
	}

	void Clear()
	{
		coefficients.clear();
	}

	static size_t SizeOfNonZero()
	{
		return sizeof(Eigen::Triplet<double>);
	}

	void Add(BigNumber i, BigNumber j, double value)
	{
		if (abs(value) > 1e-15)
			this->coefficients.push_back(Eigen::Triplet<double, SparseMatrixIndex>(i, j, value));
	}

	void Add(const NonZeroCoefficients &chunk)
	{
		this->coefficients.insert(this->coefficients.end(), chunk.coefficients.begin(), chunk.coefficients.end());
	}

	void Add(BigNumber iStart, BigNumber jStart, const DenseMatrix &m)
	{
		for (int i = 0; i < m.rows(); ++i)
		{
			for (int j = 0; j < m.cols(); ++j)
				Add(iStart + i, jStart + j, m(i, j));
		}
	}

	void Add(BigNumber iStart, BigNumber jStart, const RowMajorSparseMatrix &m)
	{
		for (int k = 0; k < m.outerSize(); ++k)
		{
			for (RowMajorSparseMatrix::InnerIterator it(m, k); it; ++it)
				Add(iStart + it.row(), jStart + it.col(), it.value());
		}
	}

	void Add(BigNumber iStart, BigNumber jStart, const ColMajorSparseMatrix &m)
	{
		for (int k = 0; k < m.outerSize(); ++k)
		{
			for (ColMajorSparseMatrix::InnerIterator it(m, k); it; ++it)
				Add(iStart + it.row(), jStart + it.col(), it.value());
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
		for (Eigen::Triplet<double, SparseMatrixIndex> t : nnzCoeffs.coefficients)
			os << t.row() << "\t" << t.col() << "\t" << t.value() << endl;
		return os;
	}
};
