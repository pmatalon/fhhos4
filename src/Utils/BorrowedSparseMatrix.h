#pragma once
#include <cassert>
#include "Types.h"

// A SparseMatrix (row-major, compressed) whose arrays belong to another code, read in place: the inputs of the library
// fhhos4 (the arrays of its Eigen matrices), which the solvers read as const SparseMatrix&. Eigen's SparseMatrix always owns its arrays:
// while it borrows the caller's, this class points the protected members of Eigen's SparseMatrix to them, and keeps its
// own (empty) arrays aside, which it gives back before Eigen's destructor runs: Eigen frees its own arrays only.
// To use as a const SparseMatrix& only: changing the matrix (resize, assignment, insertion, makeCompressed...) would
// free or reallocate the caller's arrays. Relies on the members of Eigen's SparseMatrix (m_outerSize, m_innerSize,
// m_outerIndex, m_innerNonZeros, m_data) and CompressedStorage (m_values, m_indices, m_size, m_allocatedSize), the same
// in Eigen 3.4 and 5.0 (tests/BorrowedSparseMatrixTest.cpp).
class BorrowedSparseMatrix : public SparseMatrix
{
private:
	// The non-zeros given to the matrix by Borrow(), then its own ones while it borrows the caller's
	struct Storage : public Eigen::internal::CompressedStorage<double, SparseMatrixIndex>
	{
		void Set(double* values, SparseMatrixIndex* indices, Eigen::Index size)
		{
			m_values = values;
			m_indices = indices;
			m_size = size;
			m_allocatedSize = size;
		}
	};

	Storage _aside;
	SparseMatrixIndex* _ownOuterIndex = nullptr;
	Eigen::Index _ownOuterSize = 0;
	Eigen::Index _ownInnerSize = 0;
	bool _borrowing = false;

public:
	BorrowedSparseMatrix() = default;

	// rowPointers: rows + 1 offsets in columnIndices and values, starting at 0; the columns sorted in each row
	BorrowedSparseMatrix(Eigen::Index rows, Eigen::Index cols, const SparseMatrixIndex* rowPointers, const SparseMatrixIndex* columnIndices, const double* values)
	{
		Borrow(rows, cols, rowPointers, columnIndices, values);
	}

	BorrowedSparseMatrix(const BorrowedSparseMatrix&) = delete;
	BorrowedSparseMatrix& operator=(const BorrowedSparseMatrix&) = delete;

	~BorrowedSparseMatrix()
	{
		GiveBack();
	}

	// The matrix becomes the given arrays, read in place (the previous ones given back)
	void Borrow(Eigen::Index rows, Eigen::Index cols, const SparseMatrixIndex* rowPointers, const SparseMatrixIndex* columnIndices, const double* values)
	{
		GiveBack();
		assert(!m_innerNonZeros && rowPointers && rowPointers[0] == 0);
		_ownOuterIndex = m_outerIndex;
		_ownOuterSize = m_outerSize;
		_ownInnerSize = m_innerSize;
		m_outerIndex = const_cast<SparseMatrixIndex*>(rowPointers);
		m_outerSize = rows;
		m_innerSize = cols;
		_aside.Set(const_cast<double*>(values), const_cast<SparseMatrixIndex*>(columnIndices), rowPointers[rows]);
		m_data.swap(_aside);
		_borrowing = true;
	}

	// The matrix gets its own arrays back: an empty matrix
	void GiveBack()
	{
		if (!_borrowing)
			return;
		m_data.swap(_aside);
		_aside.Set(nullptr, nullptr, 0);
		m_outerIndex = _ownOuterIndex;
		m_outerSize = _ownOuterSize;
		m_innerSize = _ownInnerSize;
		_borrowing = false;
	}
};
