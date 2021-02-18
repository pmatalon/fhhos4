#pragma once
#include "Types.h"
using namespace std;

class NNZ
{
public:
	// A must be symmetric
	static size_t StrictTriPart(const SparseMatrix& A)
	{
		return (A.nonZeros() - A.rows()) / 2;
	}
	static size_t TriPart(const SparseMatrix& A)
	{
		return StrictTriPart(A) + A.rows();
	}
};

class Cost
{
public:
	// Sparse Matrix Vector multiplication
	static size_t SpMatVec(size_t nnz)
	{
		return 2 * nnz;
	}
	// Sparse Matrix Vector multiplication
	static size_t MatVec(const SparseMatrix& A)
	{
		return SpMatVec(A.nonZeros());
	}

	// Dense dot product
	static size_t Dot(size_t rows)
	{
		return 2 * rows;
	}
	// Dense dot product
	static size_t Dot(const Vector& v)
	{
		return Dot(v.rows());
	}

	static size_t DAXPY(size_t nnz, size_t rows)
	{
		return 2 * nnz + rows;
	}
	static size_t DAXPY(const SparseMatrix& A)
	{
		return DAXPY(A.nonZeros(), A.rows());
	}

	static size_t AddVec(const Vector& v)
	{
		return v.rows();
	}

	// Sparse forward elimination
	static size_t SpFWElimination(size_t nnzL)
	{
		return 2 * nnzL;
	}

	// Sparse backward substitution
	static size_t SpBWSubstitution(size_t nnzU)
	{
		return 2 * nnzU;
	}
};