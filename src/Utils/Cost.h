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

	// AX+Y
	static size_t DAXPY(size_t nnz)
	{
		return 2 * nnz;
	}
	static size_t DAXPY(const SparseMatrix& A)
	{
		return DAXPY(A.nonZeros());
	}

	// UX+Y or LX+Y where U=upper(A), L=lower(A)
	// Requires A symmetric
	static size_t DAXPY_StrictTri(const SparseMatrix& A)
	{
		return DAXPY(NNZ::StrictTriPart(A));
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

	// Norm
	static size_t Norm(const Vector& v)
	{
		return Dot(v);
	}

	// Vector aX+Y
	static size_t VectorDAXPY(const Vector& v)
	{
		return 2 * v.rows();
	}

	// Vector x+y
	static size_t AddVec(const Vector& v)
	{
		return v.rows();
	}

	// Vector alpha*v
	static size_t ConstantVec(const Vector& v)
	{
		return v.rows();
	}


	// Sparse forward elimination (A must be symmetric)
	static size_t SpFWElimination(const SparseMatrix& A)
	{
		return SpFWElimination(NNZ::StrictTriPart(A), A.rows());
	}
	// Sparse forward elimination
	static size_t SpFWElimination(size_t nnzStrictL, size_t n)
	{
		return 2 * nnzStrictL + n;
	}

	// Sparse backward substitution (A must be symmetric)
	static size_t SpBWSubstitution(const SparseMatrix& A)
	{
		return SpBWSubstitution(NNZ::StrictTriPart(A), A.rows());
	}
	// Sparse backward substitution
	static size_t SpBWSubstitution(size_t nnzStrictU, size_t n)
	{
		return 2 * nnzStrictU + n;
	}
};