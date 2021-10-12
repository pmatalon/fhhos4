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
	static Flops SpMatVec(size_t nnz)
	{
		return 2 * nnz;
	}
	// Sparse Matrix Vector multiplication
	static Flops MatVec(const SparseMatrix& A)
	{
		return SpMatVec(A.nonZeros());
	}

	// AX+Y
	static Flops DAXPY(size_t nnz)
	{
		return 2 * nnz;
	}
	static Flops DAXPY(const SparseMatrix& A)
	{
		return DAXPY(A.nonZeros());
	}

	// UX+Y or LX+Y where U=upper(A), L=lower(A)
	// Requires A symmetric
	static Flops DAXPY_StrictTri(const SparseMatrix& A)
	{
		return DAXPY(NNZ::StrictTriPart(A));
	}

	// Dense dot product
	static Flops Dot(size_t rows)
	{
		return 2 * rows;
	}
	// Dense dot product
	static Flops Dot(const Vector& v)
	{
		return Dot(v.rows());
	}

	// Norm
	static Flops Norm(const Vector& v)
	{
		return Dot(v);
	}

	// Vector aX+Y
	static Flops VectorDAXPY(const Vector& v)
	{
		return 2 * v.rows();
	}

	// Vector x+y
	static Flops AddVec(const Vector& v)
	{
		return v.rows();
	}

	// Vector alpha*v
	static Flops ConstantVec(const Vector& v)
	{
		return v.rows();
	}


	// Sparse forward elimination (A must be symmetric)
	static Flops SpFWElimination(const SparseMatrix& A)
	{
		return SpFWElimination(NNZ::StrictTriPart(A), A.rows());
	}
	// Sparse forward elimination
	static Flops SpFWElimination(size_t nnzStrictL, size_t n)
	{
		return 2 * nnzStrictL + n;
	}

	// Sparse backward substitution (A must be symmetric)
	static Flops SpBWSubstitution(const SparseMatrix& A)
	{
		return SpBWSubstitution(NNZ::StrictTriPart(A), A.rows());
	}
	// Sparse backward substitution
	static Flops SpBWSubstitution(size_t nnzStrictU, size_t n)
	{
		return 2 * nnzStrictU + n;
	}

	static Flops LUFactorization(const SparseMatrix& A)
	{
		return 2.0 / 3.0 * pow(A.rows(), 3);
	}
	static Flops LUSolve(const SparseMatrix& L, const SparseMatrix& U)
	{
		return SpFWElimination(L) + SpBWSubstitution(U);
	}

	static Flops CholeskyFactorization(const SparseMatrix& A)
	{
		return 1.0 / 3.0 * pow(A.rows(), 3);
	}
	static Flops CholeskySolve(const SparseMatrix& L)
	{
		return SpFWElimination(L) + SpBWSubstitution(L);
	}
};