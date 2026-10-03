#pragma once
// Temporary bit-for-bit comparison of sparse matrices, to check that a new implementation gives exactly the same
// results as the old one. Not part of the build: include it in the files to instrument, then revert them, e.g.
//
//     #include "../../../../scripts/perf/BitCheck.h"   // path relative to the instrumented file
//     SparseMatrix C = SparseMatrixOps::Multiply(A, B);
//     BitCheck(C, SparseMatrix(A * B), "A*B");       // compares with the old computation
//
// Each mismatch is printed when it occurs, and the number of checks and mismatches at exit
// (grep the logs for "BITCHECK").
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include "Utils/Types.h"

inline int BitCheckCount = 0;
inline int BitCheckMismatches = 0;

// Same size, same structure (outer and inner indices) and same values, bit for bit
inline bool BitIdentical(const SparseMatrix& X0, const SparseMatrix& Y0)
{
	SparseMatrix X = X0, Y = Y0;
	X.makeCompressed();
	Y.makeCompressed();
	if (X.rows() != Y.rows() || X.cols() != Y.cols() || X.nonZeros() != Y.nonZeros())
		return false;
	return std::equal(X.outerIndexPtr(), X.outerIndexPtr() + X.outerSize() + 1, Y.outerIndexPtr())
		&& std::equal(X.innerIndexPtr(), X.innerIndexPtr() + X.nonZeros(), Y.innerIndexPtr())
		&& std::memcmp(X.valuePtr(), Y.valuePtr(), X.nonZeros() * sizeof(double)) == 0;
}

inline void BitCheck(const SparseMatrix& newResult, const SparseMatrix& oldResult, const std::string& what)
{
	BitCheckCount++;
	if (!BitIdentical(newResult, oldResult))
	{
		BitCheckMismatches++;
		std::cout << "BITCHECK MISMATCH " << what << ": nnz " << newResult.nonZeros() << " vs " << oldResult.nonZeros();
		if (newResult.rows() == oldResult.rows() && newResult.cols() == oldResult.cols())
			std::cout << ", |new - old| / |old| = " << SparseMatrix(newResult - oldResult).norm() / oldResult.norm();
		std::cout << std::endl;
	}
}

inline int _bitCheckSummaryAtExit = (std::atexit([] {
	std::cout << "BITCHECK: " << BitCheckCount << " checks, " << BitCheckMismatches << " mismatches" << std::endl;
}), 0);
