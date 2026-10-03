#pragma once
#include "../../Utils/Types.h"
#include "../../Utils/Utils.h"
using namespace std;

// Positions of the coefficients of the diagonal block in each row of a block matrix, in the arrays of the matrix
// (valuePtr(), innerIndexPtr()): in row i, [Begin[i], End[i]). The rows being sorted by column index, the coefficients
// of the diagonal block are contiguous, between those of the lower blocks and those of the upper blocks.
// Lets the smoothers skip the diagonal block without testing the column of each coefficient. With blocks of size 1,
// Begin[i] is the position of the diagonal coefficient (End[i] = Begin[i] + 1 if it is stored).
struct DiagonalBlockPositions
{
	static_assert(SparseMatrix::IsRowMajor, "The rows of the matrix must be stored contiguously.");

	vector<SparseMatrixIndex> Begin;
	vector<SparseMatrixIndex> End;

	void Setup(const SparseMatrix& A, int blockSize)
	{
		// Uncompressed, the row i doesn't end at outerIndexPtr()[i + 1]
		if (!A.isCompressed())
			Utils::FatalError("The Gauss-Seidel/SOR smoothers require a compressed matrix (see Eigen's makeCompressed()).");

		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();

		// Eigen's matrices have sorted rows (setFromTriplets(), products...), but a matrix built from external CSR arrays
		// might not: without this check, the smoothers would silently compute wrong results
		bool sortedRows = true;
		#pragma omp parallel for reduction(&&:sortedRows)
		for (BigNumber i = 0; i < A.rows(); i++)
		{
			for (SparseMatrixIndex p = outer[i] + 1; p < outer[i + 1]; p++)
				sortedRows = sortedRows && col[p - 1] < col[p];
		}
		if (!sortedRows)
			Utils::FatalError("The Gauss-Seidel/SOR smoothers require a matrix whose rows are sorted by column index, without duplicates.");

		Begin.resize(A.rows());
		End.resize(A.rows());
		#pragma omp parallel for
		for (BigNumber i = 0; i < A.rows(); i++)
		{
			SparseMatrixIndex firstCol = (i / blockSize) * blockSize;
			SparseMatrixIndex p = outer[i];
			while (p < outer[i + 1] && col[p] < firstCol)
				p++;
			Begin[i] = p;
			while (p < outer[i + 1] && col[p] < firstCol + blockSize)
				p++;
			End[i] = p;
		}
	}
};
